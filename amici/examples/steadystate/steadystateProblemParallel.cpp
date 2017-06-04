#include "steadystateProblemParallel.h"
#include "wrapfunctions.h"
#include <cstring>
#include <mpi.h>
#include <loadBalancerMaster.h>
#include <pthread.h>
#include <unistd.h>


// TODO inherit from serial
SteadystateProblemParallel::SteadystateProblemParallel(int numConditions) : numConditions(numConditions)
{
    setupUserData();
    setupExpData();

    numOptimizationParameters = udata->np;

    initialParameters = new double [numOptimizationParameters];
    fillArray(initialParameters, udata->np, 1);

    parametersMin = new double [numOptimizationParameters];
    fillArray(parametersMin, udata->np, -5);

    parametersMax = new double [numOptimizationParameters];
    fillArray(parametersMax, udata->np, 5);

    optimizationOptions = new OptimizationOptions();

    optimizationOptions->optimizer = OPTIMIZER_IPOPT;

    optimizationOptions->printToStdout = true;

    optimizationOptions->maxOptimizerIterations = 3;

    MPI_Comm_size(MPI_COMM_WORLD, &commSize);

    // generate different fixed parameter vectors
    fixedParameters = new double[udata->nk * numConditions];
    for(int i = 0; i < numConditions; ++i)
        for(int ik = 0; ik < udata->nk; ++ik)
            fixedParameters[ik + i * udata->nk] = udata->k[ik] + numConditions / 10000.0;
}

int SteadystateProblemParallel::evaluateObjectiveFunction(const double *parameters, double *objFunVal, double *objFunGrad)
{
    if(commSize > 1) {
        return evaluateParallel(parameters, objFunVal, objFunGrad);
    } else {
        return evaluateSerial(parameters, objFunVal, objFunGrad);
    }
}

int SteadystateProblemParallel::evaluateParallel(const double *parameters, double *objFunVal, double *objFunGrad)
{
    // TODO: always computes gradient; ignores simulation status

    // create load balancer job for each simulation
    JobData jobdata[numConditions];

    // mutex to wait for simulations to finish
    pthread_cond_t simulationsCond = PTHREAD_COND_INITIALIZER;
    pthread_mutex_t simulationsMutex = PTHREAD_MUTEX_INITIALIZER;

    int numJobsFinished = 0;

    for(int i = 0; i < numConditions; ++i) {
        JobData *job = &jobdata[i];
        job->jobDone = &numJobsFinished;
        job->jobDoneChangedCondition = &simulationsCond;
        job->jobDoneChangedMutex = &simulationsMutex;
        job->lenSendBuffer = sizeof(double) * (udata->nk + udata->np);
        job->sendBuffer = (char *) malloc(job->lenSendBuffer);

        double *doubleBuffer = (double *) job->sendBuffer;
        for(int ik = 0; ik < udata->nk; ++ik)
            doubleBuffer[i] = fixedParameters[udata->nk * i + ik];
        for(int ip = 0; ip < udata->np; ++ip)
            doubleBuffer[ip + udata->nk] = udata->p[ip];

        loadBalancerQueueJob(job);
    }

    // wait for simulations to finish
    pthread_mutex_lock(&simulationsMutex);
    while(numJobsFinished < numConditions) // TODO handle finished simulations here, don't wait for all to complete; stop early if errors occured
        pthread_cond_wait(&simulationsCond, &simulationsMutex);
    pthread_mutex_unlock(&simulationsMutex);
    pthread_mutex_destroy(&simulationsMutex);
    pthread_cond_destroy(&simulationsCond);

    // aggregate likelihood
    *objFunVal = 0;
    if(objFunGrad)
        fillArray(objFunGrad, numOptimizationParameters, 0.0);

    for(int i = 0; i < numConditions; ++i) {
        double *buffer = (double*) (jobdata[i].recvBuffer);
        *objFunVal -= buffer[0];

        if(objFunGrad)
            for(int ip = 0; ip < udata->np; ++ip)
                objFunGrad[ip] -= buffer[1 + ip];
        free(buffer);
    }

    sleep(1e-10);

    return 0;
}

int SteadystateProblemParallel::evaluateSerial(const double *parameters, double *objFunVal, double *objFunGrad)
{
    int status = -1;
    memcpy(udata->p, parameters, udata->np * sizeof(double));

//    printArray(parameters, udata->np);printf("\n");

    if(objFunGrad) {
        udata->sensi = AMI_SENSI_ORDER_FIRST;
        udata->sensi_meth = AMI_SENSI_FSA;
    } else {
        udata->sensi = AMI_SENSI_ORDER_NONE;
        udata->sensi_meth = AMI_SENSI_NONE;
    }

    *objFunVal = 0;
    fillArray(objFunGrad, numOptimizationParameters, 0.0);

    for(int i = 0; i < numConditions; ++i) {
        memcpy(udata->k, &fixedParameters[i * udata->nk], udata->nk * sizeof(double));

        ReturnData *rdata = getSimulationResults(udata, edata);
        status += (int) *rdata->status;

        *objFunVal -= *rdata->llh;

        if(objFunGrad)
            for(int ip = 0; ip < udata->np; ++ip)
                objFunGrad[ip] -= rdata->sllh[ip];

        delete rdata;
    }
    return status;

}

int SteadystateProblemParallel::intermediateFunction(int alg_mod, int iter_count, double obj_value, double inf_pr, double inf_du, double mu, double d_norm, double regularization_size, double alpha_du, double alpha_pr, int ls_trials)
{
    return 0;
}

void SteadystateProblemParallel::logObjectiveFunctionEvaluation(const double *parameters, double objectiveFunctionValue, const double *objectiveFunctionGradient, int numFunctionCalls, double timeElapsed)
{

}

void SteadystateProblemParallel::logOptimizerFinished(double optimalCost, const double *optimalParameters, double masterTime, int exitStatus)
{
    printf("Optimal parameters:\n\t");
    printArray(optimalParameters, udata->np);
    printf("\n");
    printf("Minimal cost: %f\n", optimalCost);
}

SteadystateProblemParallel::~SteadystateProblemParallel(){
    delete[] initialParameters;
    delete[] parametersMin;
    delete[] parametersMax;
    delete udata;
    freeExpData(edata);

    delete optimizationOptions;
}

void SteadystateProblemParallel::setupUserData()
{
    udata = new UserData(getUserData());

    udata->nt = 1;
    udata->ts = new double[udata->nt];
    udata->ts[0] = 100;

    udata->idlist = new double[udata->nx];
    fillArray(udata->idlist, udata->nx, 1);
    udata->qpositivex = new double[udata->nx];
    fillArray(udata->qpositivex, udata->nx, 1);

    udata->plist = new int[udata->np];
    udata->nplist = udata->np;
    for(int i = 0; i < udata->np; ++i) udata->plist[i] = i;

    udata->p = new double[udata->np];

    udata->k = new double[udata->nk];
    udata->k[0] = 0.1;
    udata->k[1] = 0.4;
    udata->k[2] = 0.7;
    udata->k[3] = 1;

    udata->sensi = AMI_SENSI_ORDER_FIRST;
    udata->sensi_meth = AMI_SENSI_FSA;

}

void SteadystateProblemParallel::setupExpData()
{
    edata = new ExpData();

    edata->am_my = new double[udata->nytrue * udata->nt];
    fillArray(edata->am_my, udata->nytrue * udata->nt, 1);

    edata->am_ysigma = new double[udata->nytrue * udata->nt];
    fillArray(edata->am_ysigma, udata->nytrue * udata->nt, 1);
}

