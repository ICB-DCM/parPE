#include "steadystateProblemParallel.h"
#include "wrapfunctions.h"
#include <LoadBalancerMaster.h>
#include <cstring>
#include <mpi.h>
#include <pthread.h>
#include <unistd.h>

SteadystateProblemParallel::SteadystateProblemParallel(
    parpe::LoadBalancerMaster *loadBalancer, std::string const& dataFileName)
    : loadBalancer(loadBalancer), model(std::unique_ptr<Model>(getModel())) {

    setNumOptimizationParameters(model->np);
    fillArray(initialParameters_.data(), model->np, 0);
    fillArray(parametersMin_.data(), model->np, -5);
    fillArray(parametersMax_.data(), model->np, 5);

    dataProvider = std::make_unique<SteadyStateMultiConditionDataProvider>(model.get(), dataFileName);
    udata = dataProvider->getUserData();
    numConditions = dataProvider->getNumberOfConditions();

    optimizationOptions.optimizer = parpe::OPTIMIZER_IPOPT;
    optimizationOptions.printToStdout = true;
    optimizationOptions.maxOptimizerIterations = 10;
}

int SteadystateProblemParallel::evaluateObjectiveFunction(const double *parameters, double *objFunVal, double *objFunGrad) {
    if (parpe::getMpiCommSize() > 1) {
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
    parpe::JobData jobdata[numConditions];

    // mutex to wait for simulations to finish
    pthread_cond_t simulationsCond = PTHREAD_COND_INITIALIZER;
    pthread_mutex_t simulationsMutex = PTHREAD_MUTEX_INITIALIZER;

    int numJobsFinished = 0;

    for (int i = 0; i < numConditions; ++i) {
        parpe::JobData *job = &jobdata[i];
        job->jobDone = &numJobsFinished;
        job->jobDoneChangedCondition = &simulationsCond;
        job->jobDoneChangedMutex = &simulationsMutex;
        job->sendBuffer.resize(sizeof(double) * model->np + 2 * sizeof(int));

        dataProvider->updateFixedSimulationParameters(i, *udata);
        int needGradient = objFunGrad ? 1 : 0;

        memcpy(job->sendBuffer.data(), &i, sizeof(int));
        memcpy(job->sendBuffer.data() + sizeof(int), &needGradient, sizeof(int));
        memcpy(job->sendBuffer.data() + 2 * sizeof(int), parameters,
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
        fillArray(objFunGrad, dataProvider->getNumOptimizationParameters(), 0.0);

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
        fillArray(objFunGrad, dataProvider->getNumOptimizationParameters(), 0.0);
    } else {
        udata->sensi = AMICI_SENSI_ORDER_NONE;
        udata->sensi_meth = AMICI_SENSI_NONE;
    }

    for (int i = 0; i < numConditions; ++i) {
        dataProvider->updateFixedSimulationParameters(i, *udata);
        auto edata = dataProvider->getExperimentalDataForCondition(i, udata.get());

        ReturnData *rdata = getSimulationResults(model.get(), udata.get(), edata.get());
        status += (int)*rdata->status;

        *objFunVal -= *rdata->llh;

        if (objFunGrad)
            for (int ip = 0; ip < model->np; ++ip)
                objFunGrad[ip] -= rdata->sllh[ip];

        delete rdata;
    }
    return status;
}

void SteadystateProblemParallel::messageHandler(std::vector<char> &buffer,
                                                int jobId) {
    //    int mpiRank;
    //    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    //    logmessage(LOGLVL_DEBUG, "Worker #%d: Job #%d received.", mpiRank,
    //    jobId);

    // unpack parameters

    int conditionIdx = *(int*)buffer.data();
    int needGradient = *(int*)(buffer.data() + sizeof(int));
    udata->setParameters(reinterpret_cast<double *>(buffer.data() + 2 * sizeof(int)));

    // read data for current conditions
    dataProvider->updateFixedSimulationParameters(conditionIdx, *udata);
    auto edata = dataProvider->getExperimentalDataForCondition(conditionIdx, udata.get());

    if (needGradient) {
        udata->sensi = AMICI_SENSI_ORDER_FIRST;
        udata->sensi_meth = AMICI_SENSI_FSA;
    } else {
        udata->sensi = AMICI_SENSI_ORDER_NONE;
        udata->sensi_meth = AMICI_SENSI_NONE;
    }


    // run simulation
    ReturnData *rdata = getSimulationResults(model.get(), udata.get(), edata.get());
    // printf("Result for %d: %f\n", conditionIdx, *rdata->llh);
    // pack results
    buffer.resize(sizeof(double) * (udata->nplist + 1));
    double *doubleBuffer = (double *) buffer.data();

    doubleBuffer[0] = rdata->llh[0];
    if (needGradient)
        for (int i = 0; i < udata->nplist; ++i)
            doubleBuffer[1 + i] = rdata->sllh[i];

    delete rdata;
}

SteadystateProblemParallel::~SteadystateProblemParallel() {}
