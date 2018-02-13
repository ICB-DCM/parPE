#include "standaloneSimulator.h"
#include <SimulationRunner.h>
#include <simulationResultWriter.h>
#include <amici_interface_cpp.h>
#include <LoadBalancerMaster.h>
#include <LoadBalancerWorker.h>
#include <optimizationOptions.h>

#include <iostream>

namespace parpe {

StandaloneSimulator::StandaloneSimulator(MultiConditionDataProvider *dp)
    : dataProvider(dp)
{

}


int StandaloneSimulator::run(const std::string& resultFile, const std::string& resultPath,
                             std::vector<double> parameters, LoadBalancerMaster *loadBalancer)
{
    auto edata = dataProvider->getExperimentalDataForCondition(0);
    int errors = 0;
    JobIdentifier path;

    SimulationResultWriter rw(resultFile, resultPath);
    rw.saveYMes = true;
    rw.saveYSim = true;
    rw.saveLlh = true;

    auto model = dataProvider->getModel();
    auto solver = dataProvider->getSolver();
    rw.createDatasets(*model, edata.get(), dataProvider->getNumberOfConditions());

    SimulationRunner simRunner(
                dataProvider->getNumberOfConditions(),
                [&](int simulationIdx) {
        auto myModel = std::unique_ptr<amici::Model>(model->clone());
        // extract parameters for simulation of current condition, instead
        // of sending whole  optimization parameter vector to worker
        dataProvider->updateConditionSpecificSimulationParameters(
                    simulationIdx, parameters.data(), *myModel);
        return std::pair<std::unique_ptr<amici::Model>,std::unique_ptr<amici::Solver>>(std::move(myModel), std::unique_ptr<amici::Solver>(solver->clone()));
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

        auto edata = dataProvider->getExperimentalDataForCondition(dataIdx);

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
    JobIdentifier path;
    auto model = dataProvider->getModelForCondition(0);
    auto solver = dataProvider->getSolver();
    JobAmiciSimulation<JobIdentifier> sim(solver.get(), model.get(), &path);
    sim.deserialize(buffer.data(), buffer.size());

#if QUEUE_WORKER_H_VERBOSE >= 2
    int mpiRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    printf("[%d] Received work. ", mpiRank);
    fflush(stdout);
#endif

    // do work
    JobResultAmiciSimulation result = runSimulation(path, jobId);

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

    buffer = amici::serializeToStdVec<JobResultAmiciSimulation>(result);
}


JobResultAmiciSimulation StandaloneSimulator::runSimulation(JobIdentifier path, int jobId)
{
    auto model = dataProvider->getModel();
    auto solver = dataProvider->getSolver();

    dataProvider->updateFixedSimulationParameters(path.idxConditions, *model);

    auto edata = dataProvider->getExperimentalDataForCondition(path.idxConditions);

    auto rdata = std::unique_ptr<amici::ReturnData>(
                amici::getSimulationResults(*model, edata.get(), *solver));

    RELEASE_ASSERT(rdata != nullptr, "");
    int status = *rdata->status;
    return JobResultAmiciSimulation(status, std::move(rdata), 0.0);
}



std::vector<double> getFinalParameters(std::string startIndex, H5::H5File &file)
{

    // read from last iteration (last column in /multistarts/$/iterCostFunParameters)
    std::string parameterPath = std::string("/multistarts/") + startIndex + "/iterCostFunParameters";
    H5::DataSet dataset = file.openDataSet(parameterPath);

    H5::DataSpace filespace = dataset.getSpace();
    int rank = filespace.getSimpleExtentNdims();
    RELEASE_ASSERT(rank == 2, "Rank mismatch");

    hsize_t dims[2];
    filespace.getSimpleExtentDims(dims, NULL);
    int numIter = dims[1];
    int numParam = dims[0];

    std::vector<double> parameters(numParam);

    parpe::hdf5Read2DDoubleHyperslab(file.getId(), parameterPath.c_str(),
                                     numParam, 1, 0, numIter - 1,
                                     parameters.data());
    return parameters;
}

std::vector<std::vector<double>> getParameterTrajectory(std::string startIndex, H5::H5File &file)
{
    std::string parameterPath = std::string("/multistarts/") + startIndex + "/iterCostFunParameters";
    H5::DataSet dataset = file.openDataSet(parameterPath);

    H5::DataSpace filespace = dataset.getSpace();
    int rank = filespace.getSimpleExtentNdims();
    RELEASE_ASSERT(rank == 2, "Rank mismatch");

    hsize_t dims[2];
    filespace.getSimpleExtentDims(dims, NULL);
    int numIter = dims[1];
    int numParam = dims[0];

    std::vector<std::vector<double>> parameters(numIter);

    for(int iter = 0; iter < numIter; ++iter) {
        parameters[iter] = std::vector<double>(numParam);
        parpe::hdf5Read2DDoubleHyperslab(file.getId(), parameterPath.c_str(),
                                         numParam, 1, 0, iter,
                                         parameters[iter].data());
    }

    return parameters;
}

int getNumStarts(H5::H5File file, std::string rootPath)  {
    auto o = parpe::OptimizationOptions::fromHDF5(file.getId());
    return o->numStarts;
}

int runFinalParameters(StandaloneSimulator &sim, std::string inFileName, std::string resultFileName, std::string resultPath, LoadBalancerMaster *loadBalancer) {

    H5::H5File file(inFileName, H5F_ACC_RDONLY);
    int errors = 0;

    for(int i = 0; i < getNumStarts(file); ++i) {
        std::cout<<"Running for start "<<i<<std::endl;
        try {
            auto parameters = parpe::getFinalParameters(std::to_string(i), file);
            std::string curResultPath = resultPath + "multistarts/" + std::to_string(i);
            errors += sim.run(resultFileName, curResultPath, parameters, loadBalancer);
        } catch (std::exception e) {
            std::cerr<<e.what()<<std::endl;
        }
    }

    return errors;
}

int runAlongTrajectory(StandaloneSimulator &sim, std::string inFileName, std::string resultFileName, std::string resultPath, LoadBalancerMaster *loadBalancer)
{
    H5::H5File file(inFileName, H5F_ACC_RDONLY);
    int errors = 0;

    for(int i = 0; i < getNumStarts(file); ++i) {
        try {

            auto parameters = getParameterTrajectory(std::to_string(i), file);

            for(int iter = 0; (unsigned) iter < parameters.size(); ++iter) {
                std::cout<<"Running for start "<<i<<" iter "<<iter<<std::endl;
                std::string curResultPath = resultPath + "multistarts/" + std::to_string(i) + "/iter/" + std::to_string(iter);

                errors += sim.run(resultFileName, curResultPath, parameters[iter], loadBalancer);
            }
        } catch (std::exception e) {
            std::cerr<<e.what()<<std::endl;
        }

    }

    return errors;
}

int runSimulationTasks(StandaloneSimulator &sim, std::string simulationMode,
                       std::string inFileName, std::string dataFilePath,
                       std::string resultFileName, std::string resultPath,
                       LoadBalancerMaster *loadBalancer) {

    if(simulationMode == "--at-optimum") {
        return parpe::runFinalParameters(sim, inFileName, resultFileName, resultPath, loadBalancer);
    } else if (simulationMode == "--along-trajectory") {
        return parpe::runAlongTrajectory(sim, inFileName, resultFileName, resultPath, loadBalancer);
    }

    return -1;
}

int runSimulator(MultiConditionDataProvider &dp, std::string simulationMode, std::string dataFileName, std::string dataFilePath, std::string resultFileName, std::string resultPath)
{
    parpe::StandaloneSimulator sim(&dp);
    int status = 0;
    int commSize = parpe::getMpiCommSize();

    if (commSize > 1) {
        if (parpe::getMpiRank() == 0) {
            parpe::LoadBalancerMaster loadBalancer;
            loadBalancer.run();
            status = runSimulationTasks(sim, simulationMode, dataFileName, dataFilePath, resultFileName, resultPath, &loadBalancer);
            loadBalancer.terminate();
            loadBalancer.sendTerminationSignalToAllWorkers();
        } else {
            parpe::LoadBalancerWorker lbw;
            lbw.run([&sim](std::vector<char> &buffer, int jobId) {
                sim.messageHandler(buffer, jobId);
            });
        }
    } else {
        status = runSimulationTasks(sim, simulationMode, dataFileName, dataFilePath, resultFileName, resultPath, nullptr);
    }

    return status;
}


} // namespace parpe
