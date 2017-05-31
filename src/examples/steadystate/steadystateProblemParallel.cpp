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

    numOptimizationParameters = udata->am_np;

    initialParameters = new double [numOptimizationParameters];
    fillArray(initialParameters, udata->am_np, 1);

    parametersMin = new double [numOptimizationParameters];
    fillArray(parametersMin, udata->am_np, -5);

    parametersMax = new double [numOptimizationParameters];
    fillArray(parametersMax, udata->am_np, 5);

    optimizationOptions = new OptimizationOptions();

    optimizationOptions->optimizer = OPTIMIZER_IPOPT;

    optimizationOptions->printToStdout = true;

    optimizationOptions->maxOptimizerIterations = 3;

    MPI_Comm_size(MPI_COMM_WORLD, &commSize);

    // generate different fixed parameter vectors
    fixedParameters = new double[udata->am_nk * numConditions];
    for(int i = 0; i < numConditions; ++i)
        for(int ik = 0; ik < udata->am_nk; ++ik)
            fixedParameters[ik + i * udata->am_nk] = udata->am_k[ik] + numConditions / 10000.0;
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
        job->lenSendBuffer = sizeof(double) * (udata->am_nk + udata->am_np);
        job->sendBuffer = (char *) malloc(job->lenSendBuffer);

        double *doubleBuffer = (double *) job->sendBuffer;
        for(int ik = 0; ik < udata->am_nk; ++ik)
            doubleBuffer[i] = fixedParameters[udata->am_nk * i + ik];
        for(int ip = 0; ip < udata->am_np; ++ip)
            doubleBuffer[ip + udata->am_nk] = udata->am_p[ip];

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
            for(int ip = 0; ip < udata->am_np; ++ip)
                objFunGrad[ip] -= buffer[1 + ip];
        free(buffer);
    }

    sleep(1e-10);

    return 0;
}

int SteadystateProblemParallel::evaluateSerial(const double *parameters, double *objFunVal, double *objFunGrad)
{
    int status = -1;
    memcpy(udata->am_p, parameters, udata->am_np * sizeof(double));

//    printArray(parameters, udata->am_np);printf("\n");

    if(objFunGrad) {
        udata->am_sensi = AMI_SENSI_ORDER_FIRST;
        udata->am_sensi_meth = AMI_SENSI_FSA;
    } else {
        udata->am_sensi = AMI_SENSI_ORDER_NONE;
        udata->am_sensi_meth = AMI_SENSI_NONE;
    }

    *objFunVal = 0;
    fillArray(objFunGrad, numOptimizationParameters, 0.0);

    for(int i = 0; i < numConditions; ++i) {
        memcpy(udata->am_k, &fixedParameters[i * udata->am_nk], udata->am_nk * sizeof(double));

        int tmpStatus;
        ReturnData *rdata = getSimulationResults(udata, edata, &tmpStatus);
        status += tmpStatus;

        *objFunVal -= *rdata->am_llhdata;

        if(objFunGrad)
            for(int ip = 0; ip < udata->am_np; ++ip)
                objFunGrad[ip] -= rdata->am_sllhdata[ip];

        freeReturnData(rdata);
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
    printArray(optimalParameters, udata->am_np);
    printf("\n");
    printf("Minimal cost: %f\n", optimalCost);
}

SteadystateProblemParallel::~SteadystateProblemParallel(){
    delete[] initialParameters;
    delete[] parametersMin;
    delete[] parametersMax;
    freeUserData(udata);
    freeExpData(edata);

    delete optimizationOptions;
}

void SteadystateProblemParallel::setupUserData()
{
    udata = getDefaultUserData();
    init_modeldims(udata);

    udata->am_atol = 1e-8;
    udata->am_rtol = 1e-8;

    udata->am_nt = 1;
    udata->am_ts = new double[udata->am_nt];
    udata->am_ts[0] = 100;

    udata->am_idlist = new double[udata->am_nx];
    fillArray(udata->am_idlist, udata->am_nx, 1);
    udata->am_qpositivex = new double[udata->am_nx];
    fillArray(udata->am_qpositivex, udata->am_nx, 1);

    udata->am_plist = new int[udata->am_np];
    udata->am_nplist = udata->am_np;
    for(int i = 0; i < udata->am_np; ++i) udata->am_plist[i] = i;

    udata->am_p = new double[udata->am_np];

    udata->am_k = new double[udata->am_nk];
    udata->am_k[0] = 0.1;
    udata->am_k[1] = 0.4;
    udata->am_k[2] = 0.7;
    udata->am_k[3] = 1;

    udata->am_lmm = 1;
    udata->am_iter = 1;
    udata->am_linsol = AMI_KLU;

    udata->am_maxsteps = 1e5;

    udata->am_sensi = AMI_SENSI_ORDER_FIRST;
    udata->am_sensi_meth = AMI_SENSI_FSA;

    processUserData(udata);
}

void SteadystateProblemParallel::setupExpData()
{
    edata = new ExpData();

    edata->am_my = new double[udata->am_nytrue * udata->am_nt];
    fillArray(edata->am_my, udata->am_nytrue * udata->am_nt, 1);

    edata->am_ysigma = new double[udata->am_nytrue * udata->am_nt];
    fillArray(edata->am_ysigma, udata->am_nytrue * udata->am_nt, 1);
}

