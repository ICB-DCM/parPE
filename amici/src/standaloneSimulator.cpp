#include "standaloneSimulator.h"
#include <SimulationRunner.h>
#include <simulationResultWriter.h>
#include <amici_interface_cpp.h>

namespace parpe {

StandaloneSimulator::StandaloneSimulator(MultiConditionDataProvider *dp)
    : dataProvider(dp)
{

}


int StandaloneSimulator::run(const std::string& resultFile, const std::string& resultPath,
                             std::vector<double> parameters, LoadBalancerMaster *loadBalancer)
{
    auto udata = dataProvider->getUserData();
    auto edata = dataProvider->getExperimentalDataForCondition(0, udata.get());
    int errors = 0;
    JobIdentifier path;

    SimulationResultWriter rw(resultFile, resultPath);
    rw.saveYMes = true;
    rw.saveYSim = true;
    rw.saveLlh = true;

    auto model = dataProvider->getModel();
    rw.createDatasets(*model, udata.get(), edata.get(), dataProvider->getNumberOfConditions());

    SimulationRunner simRunner(
                dataProvider->getNumberOfConditions(),
                [&](int simulationIdx) {
        // extract parameters for simulation of current condition, instead
        // of sending whole  optimization parameter vector to worker
        dataProvider->updateConditionSpecificSimulationParameters(
                    simulationIdx, parameters.data(), udata.get());
        return *udata;
    },
    [&](int simulationIdx) {
        path.idxConditions = simulationIdx;
        return path;
    },
    [&](JobData *job, int dataIdx) {

        JobResultAmiciSimulation result =
                amici::deserializeFromChar<JobResultAmiciSimulation>(
                    job->recvBuffer.data(), job->recvBuffer.size());
        job->recvBuffer = std::vector<char>(); // free buffer

        auto edata = dataProvider->getExperimentalDataForCondition(dataIdx, udata.get());

        rw.saveSimulationResults(edata.get(), result.rdata.get(), dataIdx);
    }, nullptr);


    if (loadBalancer && loadBalancer->isRunning()) {
        errors += simRunner.runDistributedMemory(loadBalancer);
    } else {
        errors += simRunner.runSharedMemory(
            [&](std::vector<char> &buffer, int jobId) {
                messageHandler(buffer, jobId);
            });
    }

    return errors;
}


void StandaloneSimulator::messageHandler(std::vector<char> &buffer, int jobId)
{
    // unpack
    auto udata = std::unique_ptr<amici::UserData>(dataProvider->getUserData());
    JobIdentifier path;
    JobAmiciSimulation<JobIdentifier> sim(udata.get(), &path);
    sim.deserialize(buffer.data(), buffer.size());

#if QUEUE_WORKER_H_VERBOSE >= 2
    int mpiRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    printf("[%d] Received work. ", mpiRank);
    fflush(stdout);
#endif

    // do work
    JobResultAmiciSimulation result = runSimulation(udata.get(), path, jobId);

#if QUEUE_WORKER_H_VERBOSE >= 2
    printf("[%d] Work done. ", mpiRank);
    fflush(stdout);
#endif

    // pack & cleanup
    delete[] result.rdata->J;
    result.rdata->J = nullptr;

    delete[] result.rdata->sigmay;
    result.rdata->sigmay = nullptr;

    delete[] result.rdata->ssigmay;
    result.rdata->ssigmay = nullptr;

    delete[] result.rdata->sx0;
    result.rdata->sx0 = nullptr;

    delete[] result.rdata->x;
    result.rdata->x = nullptr;

    delete[] result.rdata->x0;
    result.rdata->x0 = nullptr;

    delete[] result.rdata->xdot;
    result.rdata->xdot = nullptr;

    buffer = amici::serializeToStdVec<JobResultAmiciSimulation>(&result);
}


JobResultAmiciSimulation StandaloneSimulator::runSimulation(amici::UserData *udata, JobIdentifier path, int jobId)
{
    auto model = dataProvider->getModel();

    dataProvider->updateFixedSimulationParameters(path.idxConditions, *udata);

    auto edata = dataProvider->getExperimentalDataForCondition(path.idxConditions, udata);

    auto rdata = std::unique_ptr<amici::ReturnData>(
                amici::getSimulationResults(model, udata, edata.get()));

    RELEASE_ASSERT(rdata.get() != NULL, "");
    int status = *rdata->status;
    return JobResultAmiciSimulation(status, std::move(rdata), 0.0);
}



} // namespace parpe
