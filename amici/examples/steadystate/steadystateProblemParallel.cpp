#include "steadystateProblemParallel.h"
#include "wrapfunctions.h"
#include <logging.h>
#include <LoadBalancerMaster.h>
#include <cstring>
#include <mpi.h>
#include <pthread.h>
#include <unistd.h>


ExampleSteadystateGradientFunctionParallel::ExampleSteadystateGradientFunctionParallel(parpe::LoadBalancerMaster *loadBalancer, const std::string &dataFileName)
    :loadBalancer(loadBalancer),
      model(std::unique_ptr<amici::Model>(getModel()))
{
    dataProvider = std::make_unique<SteadyStateMultiConditionDataProvider>(
                std::unique_ptr<amici::Model>(model->clone()),
                dataFileName);
    numConditions = dataProvider->getNumberOfConditions();
}

parpe::FunctionEvaluationStatus ExampleSteadystateGradientFunctionParallel::evaluate(const double * const parameters, double &fval, double *gradient) const

{
    if (parpe::getMpiCommSize() > 1) {
        return evaluateParallel(parameters, fval, gradient) ? parpe::functionEvaluationFailure : parpe::functionEvaluationSuccess;
    } else {
        return evaluateSerial(parameters, fval, gradient) ? parpe::functionEvaluationFailure : parpe::functionEvaluationSuccess;
    }
}

int ExampleSteadystateGradientFunctionParallel::numParameters() const
{
    return model->np();
}

int ExampleSteadystateGradientFunctionParallel::evaluateParallel(const double *parameters,
                                                                 double &objFunVal,
                                                                 double *objFunGrad) const {
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
        job->sendBuffer.resize(sizeof(double) * model->np() + 2 * sizeof(int));

        dataProvider->updateFixedSimulationParameters(i, *model);
        int needGradient = objFunGrad ? 1 : 0;

        memcpy(job->sendBuffer.data(), &i, sizeof(int));
        memcpy(job->sendBuffer.data() + sizeof(int), &needGradient, sizeof(int));
        memcpy(job->sendBuffer.data() + 2 * sizeof(int), parameters,
               model->np() * sizeof(double));

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
    objFunVal = 0;
    if (objFunGrad)
        std::fill(objFunGrad, objFunGrad + dataProvider->getNumOptimizationParameters(), 0.0);

    for (int i = 0; i < numConditions; ++i) {
        double *buffer = (double *)(jobdata[i].recvBuffer.data());
        objFunVal -= buffer[0];

        if (objFunGrad)
            for (int ip = 0; ip < model->np(); ++ip)
                objFunGrad[ip] -= buffer[1 + ip];
    }

    sleep(1e-10);

    return 0;
}

int ExampleSteadystateGradientFunctionParallel::evaluateSerial(const double *parameters,
                                                               double &objFunVal,
                                                               double *objFunGrad) const {
    int status = 0;
    model->setParameters(std::vector<double>(parameters, parameters + numParameters()));
    //    printArray(parameters, udata->np);printf("\n");

    objFunVal = 0;

    if (objFunGrad) {
        solver->setSensitivityOrder(amici::AMICI_SENSI_ORDER_FIRST);
        solver->setSensitivityMethod(amici::AMICI_SENSI_FSA);
        std::fill_n(objFunGrad, dataProvider->getNumOptimizationParameters(), 0.0);
    } else {
        solver->setSensitivityOrder(amici::AMICI_SENSI_ORDER_NONE);
        solver->setSensitivityMethod(amici::AMICI_SENSI_NONE);
    }

    for (int i = 0; i < numConditions; ++i) {
        dataProvider->updateFixedSimulationParameters(i, *model);
        auto edata = dataProvider->getExperimentalDataForCondition(i);

        auto rdata = amici::getSimulationResults(*model, edata.get(), *solver);
        status += (int)*rdata->status;

        objFunVal -= *rdata->llh;

        if (objFunGrad)
            for (int ip = 0; ip < model->np(); ++ip)
                objFunGrad[ip] -= rdata->sllh[ip];

        delete rdata;
    }

    return status;
}

void ExampleSteadystateGradientFunctionParallel::messageHandler(std::vector<char> &buffer,
                                                                int jobId) {
    //    int mpiRank;
    //    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    //    parpe::logmessage(parpe::LOGLVL_DEBUG, "Worker #%d: Job #%d received.", mpiRank, jobId);

    // unpack parameters

    int conditionIdx = *(int*)buffer.data();
    int needGradient = *(int*)(buffer.data() + sizeof(int));
    double *pstart = reinterpret_cast<double *>(buffer.data() + 2 * sizeof(int));
    model->setParameters(std::vector<double>(pstart, pstart + model->np()));

    // read data for current conditions
    dataProvider->updateFixedSimulationParameters(conditionIdx, *model);
    auto edata = dataProvider->getExperimentalDataForCondition(conditionIdx);

    if (needGradient) {
        solver->setSensitivityOrder(amici::AMICI_SENSI_ORDER_FIRST);
        solver->setSensitivityMethod(amici::AMICI_SENSI_FSA);
    } else {
        solver->setSensitivityOrder(amici::AMICI_SENSI_ORDER_NONE);
        solver->setSensitivityMethod(amici::AMICI_SENSI_NONE);
    }

    // run simulation
    auto rdata = amici::getSimulationResults(*model, edata.get(), *solver);
    // printf("Result for %d: %f\n", conditionIdx, *rdata->llh);
    // pack results
    buffer.resize(sizeof(double) * (model->nplist() + 1));
    double *doubleBuffer = (double *) buffer.data();

    doubleBuffer[0] = rdata->llh[0];
    if (needGradient)
        for (int i = 0; i < model->nplist(); ++i)
            doubleBuffer[1 + i] = rdata->sllh[i];

    delete rdata;
}
