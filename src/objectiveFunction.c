#include "objectiveFunction.h"
#include "dataprovider.h"

#include <stdio.h>
#include <time.h>
#include <mpi.h>
#include <alloca.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>

#ifdef USE_MPE
#include <mpe.h>
#endif

#include "queueworker.h"
#include "queuemaster.h"
#include "simulationworker.h"
#include "misc.h"


#include <include/amici.h>
#include <include/udata_accessors.h>
#include <include/rdata_accessors.h>

#define MIN(a,b) (((a)<(b))?(a):(b))

#ifdef USE_MPE
// MPE event IDs for logging
extern const int mpe_event_begin_aggregate, mpe_event_end_aggregate;
extern const int mpe_event_begin_getrefs, mpe_event_end_getrefs;
extern const int mpe_event_begin_getdrugs, mpe_event_end_getdrugs;
#endif

/******************************/

// static function prototypes
static int simulateReferenceExperiments(datapath datapath, int numGenotypes, double llhRef[], double sllhRef[], UserData *udata);

static int simulateDrugExperiments(datapath path, int numGenotypes, double *llhDrug[], double *sllhDrug[], UserData *udata);

static void saveSimulationResults(datapath path);

static void updateInitialConditions(double destination[], const double src[], int count);

static bool reachedSteadyState(const double *xdot, const double *x, int numTimepoints, int numStates, double tolerance);

static void aggregateLikelihoodAndGradient(int numGenotypes, datapath path, UserData *udata, double *llhRef, double **llhDrug, double *sllhRef, double **sllhDrug, double *objectiveFunctionValue, double *objectiveFunctionGradient);

static double getLogLikelihoodIncrementForExperiment(const double llhCase, const double llhControl, const double y, const double sigmaY, double *_growthInhibition);

static void updateLogLikelihoodGradientForExperiment(const double llhCase, const double llhControl, const double *sllhCase, const double *sllhControl, const double caseY, const double caseSigmaY, const int numTheta, const double inhib, double *dloglik);

/******************************/

int evaluateObjectiveFunction(const double theta[], const int lenTheta, datapath path,
                              double *objectiveFunctionValue, double *objectiveFunctionGradient, AMI_parameter_scaling scaling) {

    UserData *udata = getMyUserData(); // TODO datapath
    udata->am_pscale = scaling;

    if(!objectiveFunctionGradient) {
        udata->am_sensi_meth = 0;
    }

    // update parameters in UserData
    assert(np == lenTheta);
    for(int i = 0; i < np; ++i) {
        p[i] = theta[i];
    }

    int numGenotypes = getNumGenotypes(path);

    // init results arrays for likelihood and likelihood sensitivity
    double llhRef[numGenotypes];
    fillArray(llhRef, numGenotypes, NAN);
    double *sllhRef = malloc(sizeof(double) * numGenotypes * np);
    fillArray(sllhRef, numGenotypes * np, NAN);

    int errors = 0;

#ifdef USE_MPE
    MPE_Log_event(mpe_event_begin_getrefs, 0, "getrefs");
#endif

    errors = simulateReferenceExperiments(path, numGenotypes, llhRef, sllhRef, udata);

#ifdef USE_MPE
    MPE_Log_event(mpe_event_end_getrefs, 0, "getrefs");
#endif

    if(errors == 0) { // only simulate drug treatments if not prior errors occured
        // result arrays
        double *llhDrug[numGenotypes];
        double *sllhDrug[numGenotypes];

        // printf("Reference simulations done, next drugs...\n"); fflush(stdout);
#ifdef USE_MPE
        MPE_Log_event(mpe_event_begin_getdrugs, 0, "getdrugs");
#endif

        errors = simulateDrugExperiments(path, numGenotypes, llhDrug, sllhDrug, udata);

#ifdef USE_MPE
        MPE_Log_event(mpe_event_end_getdrugs, 0, "getdrugs");
#endif

        // printf("Simulations done, aggregating...\n"); fflush(stdout);

        // agregate & free buffers // TODO check data file  cell line index
#ifdef USE_MPE
        MPE_Log_event(mpe_event_begin_aggregate, 0, "agg");
#endif

        aggregateLikelihoodAndGradient(numGenotypes, path, udata, llhRef, llhDrug, sllhRef, sllhDrug, objectiveFunctionValue, objectiveFunctionGradient);

#ifdef USE_MPE
        MPE_Log_event(mpe_event_end_aggregate, 0, "agg");
#endif
    }

    freeUserDataC(udata);
    free(sllhRef);

    if(errors) {
        char strBuf[100];
        sprintDatapath(strBuf, path);
        logmessage(LOGLVL_ERROR, "%s: Objective function evaluation failed!", strBuf);
    }

    return errors;
}

static void saveSimulationResults(datapath path)
{
    // TODO Implement me
    // or better save on worker side and collect later
}

static void aggregateLikelihoodAndGradient(int numGenotypes, datapath path, UserData *udata,
                                           double *llhRef, double **llhDrug, double *sllhRef, double **sllhDrug,
                                           double *objectiveFunctionValue, double *objectiveFunctionGradient)
{
    double logLikelihood = 0;
    if(objectiveFunctionGradient)
        zeros(objectiveFunctionGradient, np);

    for(int genotypeIdx = 0; genotypeIdx < numGenotypes; ++genotypeIdx) {

        path.idxGenotype = genotypeIdx + 1; // starting from 1
        const int numExperiments = getExperimentCountForCellline(path);

        for(int experimentIdx = 0; experimentIdx < numExperiments; ++experimentIdx) {

            path.idxExperiment = experimentIdx;
            ExpData *edata = getExperimentalDataForExperiment(path, udata);

            double growthInhibition;
            double logLikelihoodIncrement = getLogLikelihoodIncrementForExperiment(llhDrug[genotypeIdx][experimentIdx],
                                                                                    llhRef[genotypeIdx],
                                                                                    edata->am_my[0], edata->am_ysigma[0], &growthInhibition);
            // logmessage(LOGLVL_DEBUG, "Agg llh: old: %f inc: %f new: %f (drug: %f ref: %f) \n", logLikelihood, logLikelihoodIncrement, logLikelihood + logLikelihoodIncrement, llhDrug[genotypeIdx][experimentIdx], llhRef[genotypeIdx]);
            logLikelihood += logLikelihoodIncrement;

            if(objectiveFunctionGradient) {
                updateLogLikelihoodGradientForExperiment(llhDrug[genotypeIdx][experimentIdx], llhRef[genotypeIdx],
                                                          &sllhDrug[genotypeIdx][experimentIdx * np], &sllhRef[genotypeIdx * np],
                                                          edata->am_my[0], edata->am_ysigma[0],
                                                          np, growthInhibition, objectiveFunctionGradient);
            }
            myFreeExpData(edata);
        }

        free(llhDrug[genotypeIdx]);
        free(sllhDrug[genotypeIdx]);
    }

    // Take negative loglikelihood
    *objectiveFunctionValue = logLikelihood;
    // logmessage(LOGLVL_DEBUG, "Aggregate: llh = %e", *objectiveFunctionValue);

    if(objectiveFunctionGradient) {
        for(int i = 0; i < np; ++i) {
            objectiveFunctionGradient[i] *= 1;
            //if(fabs(objectiveFunctionGradient[i]) > 1e-18)
            //    logmessage(LOGLVL_DEBUG, "\t\t k%d %e", i, objectiveFunctionGradient[i]);
        }
    }
}


static int simulateReferenceExperiments(datapath path, int numGenotypes, double llhRef[], double sllhRef[], UserData *udata)
{
    int errors = 0;
    queueData *data = malloc(sizeof(queueData) * numGenotypes);

    int lenSendBuffer = getLengthWorkPackageMessage(np);
    int lenRecvBuffer = getLengthResultPackageMessage(np);

    char *recvBuffer = malloc(numGenotypes * lenRecvBuffer);
    char *sendBuffer = malloc(numGenotypes * lenSendBuffer);

    path.idxExperiment = EXPERIMENT_INDEX_CONTROL;

    // use recursive mutex to wait for simulations to finish
    pthread_cond_t simulationsCond = PTHREAD_COND_INITIALIZER;
    pthread_mutex_t simulationsMutex = PTHREAD_MUTEX_INITIALIZER;
    int numJobsTotal = numGenotypes;
    int numJobsFinished = 0;

    for(int celllineIdx = 0; celllineIdx < numGenotypes; ++celllineIdx) {

        path.idxGenotype = celllineIdx + 1; // starting from 1

        queueData *d = &data[celllineIdx];
        d->jobDone = &numJobsFinished;
        d->jobDoneChangedCondition = &simulationsCond;
        d->jobDoneChangedMutex = &simulationsMutex;
        d->lenRecvBuffer = lenRecvBuffer;
        d->lenSendBuffer = lenSendBuffer;
        d->sendBuffer = &sendBuffer[celllineIdx * lenSendBuffer];
        d->recvBuffer = &recvBuffer[celllineIdx * lenRecvBuffer];

        workPackageMessage work;
        work.path = path;
        work.sensitivityMethod = sensi_meth;
        work.theta = p;

        serializeWorkPackageMessage(work, np, d->sendBuffer);

        queueSimulation(d);

        // printf("Queued work: "); printDatapath(path);
    }

    // wait for simulations to finish
    pthread_mutex_lock(&simulationsMutex);
    while(numJobsFinished < numJobsTotal)
        pthread_cond_wait(&simulationsCond, &simulationsMutex);
    pthread_mutex_unlock(&simulationsMutex);
    pthread_mutex_destroy(&simulationsMutex);
    pthread_cond_destroy(&simulationsCond);

    // unpack
    for(int celllineIdx = 0; celllineIdx < numGenotypes; ++celllineIdx) {
        queueData *d = &data[celllineIdx];
        int status;

        deserializeResultPackageMessage(d->recvBuffer, np, &status, &llhRef[celllineIdx], &sllhRef[celllineIdx * np]);
        if(status != 0)
            ++errors;
    }

    free(recvBuffer);
    free(sendBuffer);
    free(data);

    return errors;
}


int simulateDrugExperiments(datapath path, int numGenotypes, double *llhDrug[], double *sllhDrug[], UserData *udata)
{
    int errors = 0;
    // TODO: need to send steadystate as initial conditions
    int lenSendBuffer = getLengthWorkPackageMessage(np);
    int lenRecvBuffer = getLengthResultPackageMessage(np);

    // count number of experiments and allocate result memory
    queueData *data[numGenotypes];
    char *recvBuffer[numGenotypes];
    char *sendBuffer[numGenotypes];

    for(int celllineIdx = 0; celllineIdx < numGenotypes; ++celllineIdx) {
        path.idxGenotype = celllineIdx + 1; // starting from 1
        int numExperiments = getExperimentCountForCellline(path);

         llhDrug[celllineIdx] = malloc(sizeof(double) * numExperiments);
        sllhDrug[celllineIdx] = malloc(sizeof(double) * numExperiments * np);

                 data[celllineIdx] = malloc(sizeof(queueData) * numExperiments);
           recvBuffer[celllineIdx] = malloc(lenRecvBuffer     * numExperiments);
           sendBuffer[celllineIdx] = malloc(lenSendBuffer     * numExperiments);
    }

    // use recursive mutex to wait for simulations to finish
    pthread_cond_t simulationsCond = PTHREAD_COND_INITIALIZER;
    pthread_mutex_t simulationsMutex = PTHREAD_MUTEX_INITIALIZER;
    int numJobsTotal = 0;
    int numJobsFinished = 0;

    // send jobs to queue
    for(int celllineIdx = 0; celllineIdx < numGenotypes; ++celllineIdx) {
        path.idxGenotype = celllineIdx + 1; // starting from 1
        int numExperiments = getExperimentCountForCellline(path);

        for(int expIdx = 0; expIdx < numExperiments; ++expIdx) {
            path.idxExperiment = expIdx;

            queueData *d = &data[celllineIdx][expIdx];
            d->jobDone = &numJobsFinished;
            d->jobDoneChangedCondition = &simulationsCond;
            d->jobDoneChangedMutex = &simulationsMutex;
            d->lenRecvBuffer = lenRecvBuffer;
            d->lenSendBuffer = lenSendBuffer;
            d->sendBuffer    = &sendBuffer[celllineIdx][expIdx * lenSendBuffer];
            d->recvBuffer    = &recvBuffer[celllineIdx][expIdx * lenRecvBuffer];

            workPackageMessage work;
            work.path = path;
            work.sensitivityMethod = sensi_meth;
            work.theta = p;
            serializeWorkPackageMessage(work, np, d->sendBuffer);

            queueSimulation(d);

            ++numJobsTotal;

            // printf("Queued work: "); printDatapath(path);
        }
    }

    // wait for simulations to finish
    pthread_mutex_lock(&simulationsMutex);
    while(numJobsFinished < numJobsTotal)
        pthread_cond_wait(&simulationsCond, &simulationsMutex); // TODO: could already queue drug simulations here for each finished controls
    pthread_mutex_unlock(&simulationsMutex);
    pthread_mutex_destroy(&simulationsMutex);
    pthread_cond_destroy(&simulationsCond);

    // unpack memory
    for(int celllineIdx = 0; celllineIdx < numGenotypes; ++celllineIdx) {
        path.idxGenotype = celllineIdx + 1; // starting from 1
        int numExperiments = getExperimentCountForCellline(path);

        for(int expIdx = 0; expIdx < numExperiments; ++expIdx) {

            queueData *d = &data[celllineIdx][expIdx];

            int status; // use!!
            deserializeResultPackageMessage(d->recvBuffer, np, &status,
                                            &llhDrug[celllineIdx][expIdx],
                                            &sllhDrug[celllineIdx][np * expIdx]);
            if(status != 0)
                ++errors;
        }

        free(data[celllineIdx]);
        free(recvBuffer[celllineIdx]);
        free(sendBuffer[celllineIdx]);
    }

    return errors;
}


bool reachedSteadyState(const double xdot[], const double x[], const int numTimepoints, const int numStates, const double tolerance) {
    for(int state = 0; state < numStates; ++state) {
        // TODO: adapt for multiple timepoints
        assert(numTimepoints == 1);
        for(int time = 0; time < numTimepoints; ++time) {
            double sensitivity = fabs(xdot[state]) / (fabs(x[state]) + tolerance);
            if(sensitivity > tolerance) {
                // logmessage(LOGLVL_DEBUG, "No steady state: %d: x %e xdot %e relxdot %e s %e\n", state,  (x[state]), (xdot[state]) , (xdot[state]) / (x[state]), sensitivity);
                return FALSE;
            }
        }
    }
    return TRUE;
}


ReturnData *getSteadystateSolutionForExperiment(datapath path, UserData *udata, int *status, ExpData **_edata, int *iterationsDone) {
    ExpData *edata = getExperimentalDataForExperiment(path, udata);

    ReturnData *rdata = getSteadystateSolution(udata, edata, status, iterationsDone);

    if(_edata) {
        *_edata = edata;
    } else {
        myFreeExpData(edata);
    }

    return rdata;
}


ReturnData *getSteadystateSolution(UserData *udata, ExpData *edata, int *status, int *iterationDone) {
    // logmessage(LOGLVL_DEBUG, "getSteadystateSolutionForExperiment");

    ReturnData *rdata;
    bool inSteadyState = FALSE;
    int iterations = 0;

    while (!inSteadyState) {
        ++iterations;

        rdata = getSimulationResults(udata, edata, status);
        // logmessage(LOGLVL_DEBUG, "llh %e", rdata->am_llhdata[0]); fflush(stdout);

        if(*status < 0) {
            error("Failed to integrate."); // TODO add dataset info, case/control, celline
            return rdata;
        }

        inSteadyState = reachedSteadyState(xdotdata, xdata, nt, nx, XDOT_REL_TOLERANCE);

        if(inSteadyState) {
            break;
        } else if(iterations >= 100) {
            logmessage(LOGLVL_WARNING, "getSteadystateSolutionForExperiment: no steady after %d iterations... aborting...", iterations);
            *status = -1;
            break;
        }

        if(iterations % 10 == 0) {
            logmessage(LOGLVL_DEBUG, "getSteadystateSolutionForExperiment: no steady state after %d iterations... trying on...", iterations);
        }

        // use previous solution as initial conditions
        updateInitialConditions(x0data, xdata, NUM_STATE_VARIABLES);

        freeReturnData(rdata);
    }
    // logmessage(LOGLVL_DEBUG, "getSteadystateSolutionForExperiment: steadystate after %d iterations", iterations); // TODO save?

    *iterationDone = iterations;

    return rdata;
}


static void updateInitialConditions(double destination[], const double src[], const int count) {
    memcpy(destination, src, count * sizeof(double));
}


double getLogLikelihoodIncrementForExperiment(const double llhCase, const double llhControl, const double y, const double sigmaY, double *_growthInhibition) {

    double growthInhibitionSimulated = llhCase / llhControl;
    double weightedError = (growthInhibitionSimulated - y) / sigmaY;

    if(_growthInhibition)
        *_growthInhibition = growthInhibitionSimulated;

    return -0.5 * weightedError * weightedError;
}


void updateLogLikelihoodGradientForExperiment(const double llhCase, const double llhControl,
                                              const double *sllhCase, const double *sllhControl,
                                              const double caseY, const double caseSigmaY,
                                              const int numTheta, const double inhib, double *dloglik)
{
    // TODO need to adapt for multiple Y

    double llhCtrlSqrd = llhControl * llhControl;
    double ySigmaSqrd = caseSigmaY * caseSigmaY;
    double weightedError = (inhib - caseY) / ySigmaSqrd;

    //logmessage(LOGLVL_DEBUG, "getDllh: llh: %e %e ^2 %e sigma2: %e inhibsim: %e inhibmes: %e", llhCase, llhControl, llhCtrlSqrd, ySigmaSqrd, inhib, caseY);

    for(int paramIdx = 0; paramIdx < numTheta; ++paramIdx) {
        double dInhib = (sllhCase[paramIdx] * llhControl - llhCase * sllhControl[paramIdx]) / llhCtrlSqrd;
        dloglik[paramIdx] -= weightedError * dInhib;

        //if(fabs(sllhCase[i]) > 1E-16 || fabs(sllhControl[i]) > 1E-16)
        //    logmessage(LOGLVL_DEBUG, "\ti = %d dInhib: %e sCase %e sCtrl %e", i,  dInhib, sllhCase[i], sllhControl[i]);
    }
}


void objectiveFunctionGradientCheck(const double theta[], const int lenTheta, datapath path, AMI_parameter_scaling scaling, const int parameterIndices[], int numParameterIndices, double epsilon)
{
    double fc = 0; // f(theta)
    double *gradient = malloc(sizeof(double) * lenTheta);

    evaluateObjectiveFunction(theta, lenTheta, path, &fc, gradient, scaling);

    double *thetaTmp = malloc(sizeof(double) * lenTheta);
    memcpy(thetaTmp, theta, sizeof(double) * lenTheta);

    printf("Index\tGradient\tfd_f\t\t(delta)\t\tfd_c\t\t(delta)\t\tfd_b\t\t(delta)\n");

    for(int i = 0; i < numParameterIndices; ++i) {
        int curInd = parameterIndices[i];
        double fb = 0, ff = 0; // f(theta + eps) , f(theta - eps)

        thetaTmp[curInd] = theta[curInd] + epsilon;
        evaluateObjectiveFunction(thetaTmp, lenTheta, path, &ff, 0, scaling);

        thetaTmp[curInd] = theta[curInd] - epsilon;
        evaluateObjectiveFunction(thetaTmp, lenTheta, path, &fb, 0, scaling);

        printf("\t\t%f\t%f\t%f\t\n", fb, fc, ff);

        double fd_f = (ff - fc) / epsilon;

        double fd_b = (fc - fb) / epsilon;

        double fd_c = (ff - fb) / (2 * epsilon);

        thetaTmp[curInd] = theta[curInd];

        printf("%d\t%f\t%f\t(%f)\t%f\t(%f)\t%f\t(%f)\n", curInd, gradient[curInd], fd_f, gradient[curInd] - fd_f, fd_c, gradient[curInd] - fd_c, fd_b, gradient[curInd] - fd_b);
    }

    free(gradient);
    free(thetaTmp);
}
