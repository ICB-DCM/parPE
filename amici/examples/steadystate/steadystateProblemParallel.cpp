#include "steadystateProblemParallel.h"
#include "wrapfunctions.h"
#include <LoadBalancerMaster.h>
#include <cstring>
#include <mpi.h>
#include <pthread.h>
#include <unistd.h>

SteadystateProblemParallel::SteadystateProblemParallel(
    LoadBalancerMaster *loadBalancer)
    : loadBalancer(loadBalancer) {
    MPI_Comm_size(MPI_COMM_WORLD, &commSize);

    numConditions = 12;
}

int SteadystateProblemParallel::evaluateObjectiveFunction(
    const double *parameters, double *objFunVal, double *objFunGrad) {
    if (commSize > 1) {
        return evaluateParallel(parameters, objFunVal, objFunGrad);
    } else {
        return evaluateSerial(parameters, objFunVal, objFunGrad);
    }
}

int SteadystateProblemParallel::evaluateParallel(const double *parameters,
                                                 double *objFunVal,
                                                 double *objFunGrad) {
    // TODO: always computes gradient; ignores simulation status

    // create load balancer job for each simulation
    JobData jobdata[numConditions];

    // mutex to wait for simulations to finish
    pthread_cond_t simulationsCond = PTHREAD_COND_INITIALIZER;
    pthread_mutex_t simulationsMutex = PTHREAD_MUTEX_INITIALIZER;

    int numJobsFinished = 0;

    for (int i = 0; i < numConditions; ++i) {
        JobData *job = &jobdata[i];
        job->jobDone = &numJobsFinished;
        job->jobDoneChangedCondition = &simulationsCond;
        job->jobDoneChangedMutex = &simulationsMutex;
        job->lenSendBuffer = sizeof(double) * model->np + 2 * sizeof(int);
        job->sendBuffer = new char[job->lenSendBuffer];

        readFixedParameters(i);
        int needGradient = objFunGrad ? 1 : 0;

        memcpy(job->sendBuffer, &i, sizeof(int));
        memcpy(job->sendBuffer + sizeof(int), &needGradient, sizeof(int));
        memcpy(job->sendBuffer + 2 * sizeof(int), parameters,
               model->np * sizeof(double));

        loadBalancer->queueJob(job);
    }

    // wait for simulations to finish
    pthread_mutex_lock(&simulationsMutex);
    while (numJobsFinished < numConditions) // TODO handle finished simulations
                                            // here, don't wait for all to
                                            // complete; stop early if errors
                                            // occured
        pthread_cond_wait(&simulationsCond, &simulationsMutex);
    pthread_mutex_unlock(&simulationsMutex);
    pthread_mutex_destroy(&simulationsMutex);
    pthread_cond_destroy(&simulationsCond);

    // aggregate likelihood
    *objFunVal = 0;
    if (objFunGrad)
        fillArray(objFunGrad, numOptimizationParameters, 0.0);

    for (int i = 0; i < numConditions; ++i) {
        double *buffer = (double *)(jobdata[i].recvBuffer);
        *objFunVal -= buffer[0];

        if (objFunGrad)
            for (int ip = 0; ip < model->np; ++ip)
                objFunGrad[ip] -= buffer[1 + ip];
        delete[] buffer;
    }

    sleep(1e-10);

    return 0;
}

int SteadystateProblemParallel::evaluateSerial(const double *parameters,
                                               double *objFunVal,
                                               double *objFunGrad) {
    int status = 0;
    udata->setParameters(parameters);
    //    printArray(parameters, udata->np);printf("\n");

    *objFunVal = 0;

    if (objFunGrad) {
        udata->sensi = AMICI_SENSI_ORDER_FIRST;
        udata->sensi_meth = AMICI_SENSI_FSA;
        fillArray(objFunGrad, numOptimizationParameters, 0.0);
    } else {
        udata->sensi = AMICI_SENSI_ORDER_NONE;
        udata->sensi_meth = AMICI_SENSI_NONE;
    }

    for (int i = 0; i < numConditions; ++i) {
        readFixedParameters(i);
        readMeasurement(i);

        ReturnData *rdata = getSimulationResults(model, udata, edata);
        status += (int)*rdata->status;

        *objFunVal -= *rdata->llh;

        if (objFunGrad)
            for (int ip = 0; ip < model->np; ++ip)
                objFunGrad[ip] -= rdata->sllh[ip];

        delete rdata;
    }
    return status;
}

void SteadystateProblemParallel::messageHandler(char **buffer, int *size,
                                                int jobId) {
    //    int mpiRank;
    //    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    //    logmessage(LOGLVL_DEBUG, "Worker #%d: Job #%d received.", mpiRank,
    //    jobId);

    // unpack parameters
    int conditionIdx = (int)**buffer;
    int needGradient = (int)*(*buffer + sizeof(int));
    udata->setParameters(reinterpret_cast<double *>(*buffer + 2 * sizeof(int)));
    delete[] * buffer;

    // read data for current conditions
    readFixedParameters(conditionIdx);
    readMeasurement(conditionIdx);
    requireSensitivities(needGradient);

    // run simulation
    ReturnData *rdata = getSimulationResults(model, udata, edata);
    // printf("Result for %d: %f\n", conditionIdx, *rdata->llh);
    // pack results
    *size = sizeof(double) * (udata->nplist + 1);
    *buffer = new char[*size];
    double *doubleBuffer = (double *)*buffer;

    doubleBuffer[0] = rdata->llh[0];
    if (needGradient)
        for (int i = 0; i < udata->nplist; ++i)
            doubleBuffer[1 + i] = rdata->sllh[i];

    delete rdata;
}

SteadystateProblemParallel::~SteadystateProblemParallel() {}
