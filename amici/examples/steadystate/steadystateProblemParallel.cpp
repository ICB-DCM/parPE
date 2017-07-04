#include "steadystateProblemParallel.h"
#include "wrapfunctions.h"
#include <cstring>
#include <mpi.h>
#include <loadBalancerMaster.h>
#include <pthread.h>
#include <unistd.h>

SteadystateProblemParallel::SteadystateProblemParallel()
{
    MPI_Comm_size(MPI_COMM_WORLD, &commSize);

    numConditions = 12;
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
        job->lenSendBuffer = sizeof(double) * udata->np + 2 * sizeof(int);
        job->sendBuffer = (char *) malloc(job->lenSendBuffer);

        readFixedParameters(i);
        int needGradient = objFunGrad ? 1 : 0;

        memcpy(job->sendBuffer, &i, sizeof(int));
        memcpy(job->sendBuffer + sizeof(int), &needGradient, sizeof(int));
        memcpy(job->sendBuffer + 2 * sizeof(int), parameters, udata->np * sizeof(double));

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
    int status = 0;
    memcpy(udata->p, parameters, udata->np * sizeof(double));

//    printArray(parameters, udata->np);printf("\n");

    *objFunVal = 0;

    if(objFunGrad) {
        udata->sensi = AMICI_SENSI_ORDER_FIRST;
        udata->sensi_meth = AMICI_SENSI_FSA;
        fillArray(objFunGrad, numOptimizationParameters, 0.0);
    } else {
        udata->sensi = AMICI_SENSI_ORDER_NONE;
        udata->sensi_meth = AMICI_SENSI_NONE;
    }


    for(int i = 0; i < numConditions; ++i) {
        readFixedParameters(i);
        readMeasurement(i);

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


SteadystateProblemParallel::~SteadystateProblemParallel(){
}

