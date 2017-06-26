#include "multiConditionProblem.h"
#include "simulationWorkerAmici.h"
#include "multiConditionProblemResultWriter.h"

#include <math.h>
#include <string.h>
#include <assert.h>

#include <loadBalancerMaster.h>
#include <loadBalancerWorker.h>
#include <misc.h>
#include <logging.h>
#include "steadystateSimulator.h"

// For debugging:
// skip objective function evaluation completely
//#define NO_OBJ_FUN_EVAL

ReturnData *getDummyRdata(UserData *udata, int *iterationsDone);


void handleWorkPackage(char **buffer, int *msgSize, int jobId, void *userData)
{    
    MultiConditionProblem *problem = (MultiConditionProblem *) userData;
    MultiConditionDataProvider *dataProvider = problem->getDataProvider();

    // unpack
    UserData *udata = dataProvider->getUserData();
    JobIdentifier path;
    JobAmiciSimulation::toUserData(*buffer, udata, &path);
    free(*buffer);

#if QUEUE_WORKER_H_VERBOSE >= 2
    int mpiRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    printf("[%d] Received work. ", mpiRank); printDatapath(path); fflush(stdout);
#endif
    //TODO: need resultwriter here

    // work
    int status = 0;
    ReturnData *rdata = MultiConditionProblem::runAndLogSimulation(udata, dataProvider, path, jobId, NULL, &status);

#if QUEUE_WORKER_H_VERBOSE >= 2
    printf("[%d] Work done. ", mpiRank); printDatapath(path); fflush(stdout);
#endif

    // pack & cleanup
    *msgSize = JobResultAmiciSimulation::getLength(udata->np);
    *buffer = (char*) malloc(*msgSize);
    JobResultAmiciSimulation::serialize(rdata, udata, status, *buffer);

    delete rdata;
    delete udata;
}


MultiConditionProblem::MultiConditionProblem() : OptimizationProblem() // for testing only
{
    udata = NULL;
    dataProvider = NULL;

    lastOptimizationParameters = NULL;
    lastObjectiveFunctionGradient = NULL;
    lastObjectiveFunctionValue = NAN;
    path = {0};
}

MultiConditionProblem::MultiConditionProblem(MultiConditionDataProvider *dataProvider) : MultiConditionProblem()
{
    this->dataProvider = dataProvider;
    udata = dataProvider->getUserDataForCondition(0);

    if(udata == NULL)
        abort();

    init();
}

void MultiConditionProblem::init()
{  
    numOptimizationParameters = dataProvider->getNumOptimizationParameters();
    parametersMin = new double[numOptimizationParameters];
    dataProvider->getOptimizationParametersLowerBounds(parametersMin);
    parametersMax = new double[numOptimizationParameters];
    dataProvider->getOptimizationParametersUpperBounds(parametersMax);
    initialParameters = new double[numOptimizationParameters];
    getRandomStartingpoint(parametersMin, parametersMax, numOptimizationParameters, initialParameters);

    lastOptimizationParameters = new double[numOptimizationParameters];
    lastObjectiveFunctionGradient = new double[numOptimizationParameters];
}

int MultiConditionProblem::evaluateObjectiveFunction(const double *optimiziationVariables,
                                                double *objectiveFunctionValue,
                                                double *objectiveFunctionGradient)
{
#ifdef NO_OBJ_FUN_EVAL
    if(objectiveFunctionGradient)
        for(int i = 0; i < numOptimizationParameters; ++i)
            objectiveFunctionGradient = 0;
    *objectiveFunctionValue = 1;
    return 0;
#endif
    updateUserData(optimiziationVariables, objectiveFunctionGradient);

    *objectiveFunctionValue = 0;

    if(objectiveFunctionGradient)
        zeros(objectiveFunctionGradient, numOptimizationParameters);


    int errors = runSimulations(optimiziationVariables, objectiveFunctionValue, objectiveFunctionGradient);

    if(errors) {
        printObjectiveFunctionFailureMessage();
        *objectiveFunctionValue = INFINITY;
    }

    lastObjectiveFunctionValue = *objectiveFunctionValue;
    memcpy(lastOptimizationParameters, optimiziationVariables, sizeof(double) * numOptimizationParameters);
    if(objectiveFunctionGradient)
        memcpy(lastObjectiveFunctionGradient, objectiveFunctionGradient, sizeof(double) * numOptimizationParameters);

    return errors;
}

int MultiConditionProblem::intermediateFunction(int alg_mod, int iter_count, double obj_value, double inf_pr, double inf_du, double mu, double d_norm, double regularization_size, double alpha_du, double alpha_pr, int ls_trials)
{
    bool stop = false;

    path.idxLocalOptimizationIteration = iter_count;

    char strBuf[50];
    sprintJobIdentifier(strBuf, path);
    //    logmessage(LOGLVL_INFO, "%s: %d %d %e %e %e %e %e %e %e %e %d", strBuf,
    //               alg_mod, iter_count, obj_value, inf_pr, inf_du,
    //               mu, d_norm, regularization_size, alpha_du, alpha_pr, ls_trials);

    if(resultWriter)
        resultWriter->logLocalOptimizerIteration(iter_count, lastOptimizationParameters, numOptimizationParameters, obj_value, lastObjectiveFunctionGradient, 0,
                               alg_mod, inf_pr, inf_du, mu, d_norm, regularization_size, alpha_du, alpha_pr, ls_trials);
    return stop;
}

void MultiConditionProblem::logObjectiveFunctionEvaluation(const double *parameters, double objectiveFunctionValue, const double *objectiveFunctionGradient, int numFunctionCalls, double timeElapsed)
{
    if(resultWriter)
        resultWriter->logLocalOptimizerObjectiveFunctionEvaluation(parameters, numOptimizationParameters, objectiveFunctionValue, objectiveFunctionGradient, numFunctionCalls, timeElapsed);
}

void MultiConditionProblem::logOptimizerFinished(double optimalCost, const double *optimalParameters, double masterTime, int exitStatus)
{
    char strBuf[100];
    sprintJobIdentifier(strBuf, path);
    logmessage(LOGLVL_INFO, "%s: Ipopt status %d, final llh: %e, time: %f.", strBuf, exitStatus, optimalCost, masterTime);

    //    saveLocalOptimizerResults(this, path, optimalCost, optimalParameters, masterTime, exitStatus);
}

ReturnData *MultiConditionProblem::runAndLogSimulation(UserData *udata, MultiConditionDataProvider *dataProvider, JobIdentifier path, int jobId, MultiConditionProblemResultWriter *resultWriter, int *status)
{
    double startTime = MPI_Wtime();

    // run simulation
    int iterationsUntilSteadystate = 0;

    ExpData *edata = dataProvider->getExperimentalDataForExperimentAndUpdateUserData(path.idxConditions, udata);

    if(edata == NULL) {
        logmessage(LOGLVL_CRITICAL, "Failed to get experiment data. Check data file. Aborting.");
        abort();
    }

    ReturnData *rdata = getSimulationResults(udata, edata);

    freeExpData(edata);

    double endTime = MPI_Wtime();
    double timeSeconds = (endTime - startTime);

    if(isnan(rdata->llh[0]) && *status == 0)
        *status = -1;

    char pathStrBuf[100];
    sprintJobIdentifier(pathStrBuf, path);
    logmessage(LOGLVL_DEBUG, "Result for %s: %e  (%d) (%.2fs)", pathStrBuf, rdata->llh[0], *status, timeSeconds);


    // check for NaNs
    if(udata->sensi >= AMI_SENSI_ORDER_FIRST)
        for(int i = 0; i < udata->np; ++i)
            if(isnan(rdata->sllh[i]))
                logmessage(LOGLVL_DEBUG, "Result for %s: contains NaN at %d", pathStrBuf, i);
            else if(isinf(rdata->sllh[i]))
                logmessage(LOGLVL_DEBUG, "Result for %s: contains Inf at %d", pathStrBuf, i);

    // assert(*status == 0);
    // TODO write simulation results (set jobid as attribute)

    // TODO save Y
    if(resultWriter)
        resultWriter->logSimulation(udata->p, rdata->llh[0], rdata->sllh,
            timeSeconds, udata->np, udata->nx, rdata->x, rdata->sx, rdata->y,
            jobId, iterationsUntilSteadystate);

    return rdata;
}


MultiConditionProblem::~MultiConditionProblem()
{
    delete udata;

    delete[] initialParameters;
    delete[] parametersMax;
    delete[] parametersMin;
    delete[] lastObjectiveFunctionGradient;
    delete[] lastOptimizationParameters;
}

void MultiConditionProblem::updateUserData(const double *simulationParameters, const double *objectiveFunctionGradient)
{
    setSensitivityOptions(objectiveFunctionGradient);
    // update common parameters in UserData, cell-line specific ones are updated later
    for(int i = 0; i < udata->np; ++i) {
        udata->p[i] = simulationParameters[i];
    }

}

int MultiConditionProblem::runSimulations(const double *optimizationVariables, double *logLikelihood, double *objectiveFunctionGradient)
{   
    JobIdentifier path = this->path;

    int numJobsTotal = dataProvider->getNumberOfConditions();
    int numJobsFinished = 0;

    // TODO: allocate and free piecewise or according to max queue length
    JobData *data = new JobData[numJobsTotal];

    // TODO: need to send previous steadystate as initial conditions
    int lenSendBuffer = JobAmiciSimulation::getLength(udata->np, sizeof(JobIdentifier));

    // mutex to wait for simulations to finish
    pthread_cond_t simulationsCond = PTHREAD_COND_INITIALIZER;
    pthread_mutex_t simulationsMutex = PTHREAD_MUTEX_INITIALIZER;

    for(int simulationIdx = 0; simulationIdx < numJobsTotal; ++simulationIdx) {
        path.idxConditions = simulationIdx;
        updateUserDataConditionSpecificParameters(simulationIdx, optimizationVariables);
        queueSimulation(path, &data[simulationIdx],  &numJobsFinished,
                        &simulationsCond, &simulationsMutex,
                        lenSendBuffer);
        // printf("Queued work: "); printDatapath(path);
    }

    // wait for simulations to finish
    pthread_mutex_lock(&simulationsMutex);
    while(numJobsFinished < numJobsTotal) // TODO handle finished simulations here, don't wait for all to complete; stop early if errors occured
        pthread_cond_wait(&simulationsCond, &simulationsMutex);
    pthread_mutex_unlock(&simulationsMutex);
    pthread_mutex_destroy(&simulationsMutex);
    pthread_cond_destroy(&simulationsCond);

    // unpack
    int errors = aggregateLikelihood(data, logLikelihood, objectiveFunctionGradient);;
    delete[] data;

    return errors;

}

int MultiConditionProblem::aggregateLikelihood(JobData *data, double *logLikelihood, double *objectiveFunctionGradient)
{
    int errors = 0;
    double llhTmp;
    double sllhTmp[udata->np];
    int numJobsTotal = dataProvider->getNumberOfConditions();

    for(int simulationIdx = 0; simulationIdx < numJobsTotal; ++simulationIdx) {
        errors += unpackSimulationResult(&data[simulationIdx], sllhTmp, &llhTmp);

        *logLikelihood -= llhTmp;

        if(objectiveFunctionGradient)
            addSimulationGradientToObjectiveFunctionGradient(simulationIdx,
                        sllhTmp, objectiveFunctionGradient);
    }

    return errors;
}

int MultiConditionProblemSerial::runSimulations(const double *optimizationVariables, double *logLikelihood, double *objectiveFunctionGradient)
{
    JobIdentifier path = this->path;

    int numJobsTotal = dataProvider->getNumberOfConditions();

    // TODO: allocate and free piecewise or according to max queue length
    JobData *data = new JobData[numJobsTotal];

    // TODO: need to send previous steadystate as initial conditions
    int lenSendBuffer = JobAmiciSimulation::getLength(udata->np, sizeof(JobIdentifier));

    for(int simulationIdx = 0; simulationIdx < numJobsTotal; ++simulationIdx) {
        path.idxConditions = simulationIdx;

        updateUserDataConditionSpecificParameters(simulationIdx, optimizationVariables);

        data[simulationIdx].lenSendBuffer = lenSendBuffer;
        data[simulationIdx].sendBuffer = (char *) malloc(lenSendBuffer); // malloc, because will be free()'d by queue

        JobAmiciSimulation work;
        work.data = &path;
        work.lenData = sizeof(path);
        work.sensitivityMethod = udata->sensi_meth;
        work.numSimulationParameters = udata->np;
        work.simulationParameters = udata->p;

        work.serialize(data[simulationIdx].sendBuffer);
        int sizeDummy;
        handleWorkPackage(&data[simulationIdx].sendBuffer, &sizeDummy, simulationIdx, this);
        data[simulationIdx].recvBuffer = data[simulationIdx].sendBuffer;

    }

    // unpack
    int errors = aggregateLikelihood(data, logLikelihood, objectiveFunctionGradient);;

    delete[] data;


    return errors;

}

void MultiConditionProblem::printObjectiveFunctionFailureMessage()
{
    char strBuf[100];
    sprintJobIdentifier(strBuf, path);
    logmessage(LOGLVL_ERROR, "%s: Objective function evaluation failed!", strBuf);

}

void MultiConditionProblem::addSimulationGradientToObjectiveFunctionGradient(int conditionIdx, const double *simulationGradient, double *objectiveFunctionGradient)
{
    int numCommon = dataProvider->getNumCommonParameters();

    // global parameters: simply add
    for(int paramIdx = 0; paramIdx < numCommon; ++paramIdx)
        objectiveFunctionGradient[paramIdx] -= simulationGradient[paramIdx];

    int numConditionSpecificParams = dataProvider->getNumConditionSpecificParametersPerSimulation();

    // cellline specific parameters: map simulation to optimization parameters
    int firstIndexOfCurrentConditionsSpecificOptimizationParameters = dataProvider->getIndexOfFirstConditionSpecificOptimizationParameter(conditionIdx);
    for(int paramIdx = 0; paramIdx < numConditionSpecificParams; ++paramIdx) {
        int idxGrad = firstIndexOfCurrentConditionsSpecificOptimizationParameters + paramIdx;
        int idxSim  = numCommon + paramIdx;
        objectiveFunctionGradient[idxGrad] -= simulationGradient[idxSim];
    }
}

int MultiConditionProblem::unpackSimulationResult(JobData *d, double *sllhBuffer, double *llh)
{
    JobResultAmiciSimulation result;
    result.sllh = sllhBuffer;

    result.deserialize(d->recvBuffer);
    *llh = result.llh;

    free(d->recvBuffer);

    return result.status != 0;
}

void MultiConditionProblem::queueSimulation(JobIdentifier path, JobData *d, int *jobDone, pthread_cond_t *jobDoneChangedCondition, pthread_mutex_t *jobDoneChangedMutex, int lenSendBuffer)
{
    *d = initJobData(lenSendBuffer, NULL, jobDone, jobDoneChangedCondition, jobDoneChangedMutex);

    JobAmiciSimulation work;
    work.data = &path;
    work.lenData = sizeof(path);
    work.sensitivityMethod = udata->sensi_meth;
    work.numSimulationParameters = udata->np;
    work.simulationParameters = udata->p;
    work.serialize(d->sendBuffer);

    loadBalancerQueueJob(d);
}

void MultiConditionProblem::updateUserDataConditionSpecificParameters(int conditionIndex, const double *optimizationParams)
{
    /* Optimization parameters are [commonParameters, condition1SpecificParameters, condition2SpecificParameters, ...]
     * number of condition specific parameters is the same for all cell lines.
     * Simulation parameters are [commonParameters, currentCelllineSpecificParameter]
     */

    const int numCommonParams = dataProvider->getNumCommonParameters();
    const int numSpecificParams = dataProvider->getNumConditionSpecificParametersPerSimulation();

    // beginning of condiftion specific simulation parameters within optimization parameters
    const double *pConditionSpecificOptimization = &optimizationParams[dataProvider->getIndexOfFirstConditionSpecificOptimizationParameter(conditionIndex)];

    // beginning of condiftion specific simulation parameters within simulation parameters
    double *pConditionSpecificSimulation = &(udata->p[numCommonParams]);

    memcpy(pConditionSpecificSimulation, pConditionSpecificOptimization, numSpecificParams * sizeof(double));
}

void MultiConditionProblem::setSensitivityOptions(bool sensiRequired)
{
    // sensitivities requested?
    if(sensiRequired) {
        udata->sensi = AMI_SENSI_ORDER_FIRST;
        udata->sensi_meth = AMI_SENSI_ASA;
    } else {
        udata->sensi = AMI_SENSI_ORDER_NONE;
        udata->sensi_meth = AMI_SENSI_NONE;
    }
}

MultiConditionDataProvider *MultiConditionProblem::getDataProvider()
{
    return dataProvider;
}

OptimizationProblem *MultiConditionProblemGeneratorForMultiStart::getLocalProblemImpl(int multiStartIndex)
{
    int mpiCommSize = 1;
    int mpiInitialized = 0;

    MPI_Initialized(&mpiInitialized);
    if(mpiInitialized)
        MPI_Comm_size(MPI_COMM_WORLD, &mpiCommSize);

    MultiConditionProblem *problem;

    if(mpiCommSize == 1)
        problem = new MultiConditionProblemSerial(dp);
    else
        problem = new MultiConditionProblem(dp);

    problem->optimizationOptions = new OptimizationOptions(*options);
    problem->resultWriter = resultWriter;
    problem->path.idxLocalOptimization = multiStartIndex;

    return problem;
}
