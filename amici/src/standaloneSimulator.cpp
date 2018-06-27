#include "standaloneSimulator.h"
#include "SimulationRunner.h"
#include "simulationResultWriter.h"
#include <LoadBalancerMaster.h>
#include <LoadBalancerWorker.h>
#include <optimizationOptions.h>
#include <misc.h>
#include "hierarchicalOptimization.h"

#include <gsl/gsl-lite.hpp>

#include <iostream>

namespace parpe {

StandaloneSimulator::StandaloneSimulator(MultiConditionDataProvider *dp)
    : dataProvider(dp)
{

}


int StandaloneSimulator::run(const std::string& resultFile,
                             const std::string& resultPath,
                             std::vector<double> const& optimizationParameters,
                             LoadBalancerMaster *loadBalancer,
                             H5::H5File const& inputFile)
{
    // std::cout<<"file: "<<resultFile<<" path: "<<resultPath<<" lbm:"<<loadBalancer<<std::endl;
    int errors = 0;
    JobIdentifier path;

    SimulationResultWriter rw(resultFile, resultPath);
    rw.saveYMes = true;
    rw.saveYSim = true;
    rw.saveLlh = true;

    auto model = dataProvider->getModel();
    auto solver = dataProvider->getSolver();
    solver->setSensitivityOrder(amici::AMICI_SENSI_ORDER_NONE);

    std::vector<double> parameters = optimizationParameters;
    HierachicalOptimizationWrapper hierarchical(nullptr, 0, 0, 0);

    auto options = OptimizationOptions::fromHDF5(inputFile.getId(), "/inputData/optimizationOptions");
    /* if hierarchical optimization is selected and the analytical parameters have not been saved,
     * we need to recompute them. otherwise we don't need to distinguish
     */
    if(parameters.size() == (unsigned)dataProvider->getNumOptimizationParameters())
        options->hierarchicalOptimization  = false;

    if(options->hierarchicalOptimization) {
        // TODO: get rid of that. we want fun.evaluate(), independently of hierarchical or not
        auto hierarchicalScalingReader = new AnalyticalParameterHdf5Reader (inputFile,
                                                                            "/inputData/scalingParameterIndices",
                                                                            "/inputData/scalingParametersMapToObservables");
        auto hierarchicalOffsetReader = new AnalyticalParameterHdf5Reader (inputFile,
                                                                           "/inputData/offsetParameterIndices",
                                                                           "/inputData/offsetParametersMapToObservables");
        auto hierarchicalSigmaReader = new AnalyticalParameterHdf5Reader (inputFile,
                                                                           "/inputData/sigmaParameterIndices",
                                                                           "/inputData/sigmaParametersMapToObservables");
        auto proportionalityFactorIndices = hierarchicalScalingReader->getOptimizationParameterIndices();
        auto offsetParameterIndices = hierarchicalOffsetReader->getOptimizationParameterIndices();
        auto sigmaParameterIndices = hierarchicalSigmaReader->getOptimizationParameterIndices();
        auto wrappedFun = std::make_unique<AmiciSummedGradientFunction<int>>(dataProvider, loadBalancer);

        hierarchical = HierachicalOptimizationWrapper(std::move(wrappedFun),
                                                std::unique_ptr<AnalyticalParameterHdf5Reader>(hierarchicalScalingReader),
                                                std::unique_ptr<AnalyticalParameterHdf5Reader>(hierarchicalOffsetReader),
                                                std::unique_ptr<AnalyticalParameterHdf5Reader>(hierarchicalSigmaReader),
                                                dataProvider->getNumberOfConditions(), model->nytrue, model->nt(),
                                                ErrorModel::normal);
    
        // expand parameter vector
        auto scalingDummy = hierarchical.getDefaultScalingFactors();
        auto offsetDummy = hierarchical.getDefaultOffsetParameters();
        auto sigmaDummy = hierarchical.getDefaultSigmaParameters();
        parameters = spliceParameters(gsl::make_span(optimizationParameters.data(), optimizationParameters.size()),
                                      proportionalityFactorIndices, offsetParameterIndices, sigmaParameterIndices,
                                      scalingDummy, offsetDummy, sigmaDummy);

        // get outputs, scale
        // TODO need to pass aggreate function for writing
    } else {
        // is already the correct length
        // parameters = optimizationParameters;
    }

    RELEASE_ASSERT(parameters.size() == (unsigned)dataProvider->getNumOptimizationParameters(), "Size of supplied parameter vector does not match model dimensions.");

    rw.createDatasets(*model, dataProvider->getNumberOfConditions());

    SimulationRunner simRunner(
                dataProvider->getNumberOfConditions(),
                [&](int simulationIdx) { /* model & solver */
        auto myModel = std::unique_ptr<amici::Model>(model->clone());
        // extract parameters for simulation of current condition, instead
        // of sending whole  optimization parameter vector to worker
        dataProvider->updateSimulationParameters(
                    simulationIdx, parameters, *myModel);
        return std::pair<std::unique_ptr<amici::Model>,std::unique_ptr<amici::Solver>>(std::move(myModel), std::unique_ptr<amici::Solver>(solver->clone()));
    },
    [&](int simulationIdx) { /* condition id */
        path.idxConditions = simulationIdx;
        return path;
    }, [&](JobData *job, int dataIdx) { /* job finished */
        if(options->hierarchicalOptimization)
            return;

        JobResultAmiciSimulation result =
                amici::deserializeFromChar<JobResultAmiciSimulation>(
                    job->recvBuffer.data(), job->recvBuffer.size());
        job->recvBuffer = std::vector<char>(); // free buffer

        auto edata = dataProvider->getExperimentalDataForCondition(dataIdx);

        rw.saveSimulationResults(edata.get(), result.rdata.get(), dataIdx);
    },
    [&](std::vector<JobData> &jobs) -> int { /* all finished */
        if(!options->hierarchicalOptimization)
            return 0; // Work was already done in above function

        // must wait for all jobs to finish because of hierarchical optimization and scaling factors
        std::vector<JobResultAmiciSimulation> jobResults(jobs.size());
        std::vector<std::vector<double> > modelOutputs(jobs.size());

        // collect all model outputs
        for(int dataIdx = 0; (unsigned) dataIdx < jobs.size(); ++dataIdx) {
            auto& job = jobs[dataIdx];
            JobResultAmiciSimulation result = amici::deserializeFromChar<JobResultAmiciSimulation>(
                        job.recvBuffer.data(), job.recvBuffer.size());
            swap(jobResults[dataIdx], result);
            job.recvBuffer = std::vector<char>(); // free buffer
            modelOutputs[dataIdx] = jobResults[dataIdx].rdata->y;
        }

        //  compute scaling factors and offset parameters
        auto allMeasurements = dataProvider->getAllMeasurements();

        if(options->hierarchicalOptimization) {
            auto scalings = hierarchical.computeAnalyticalScalings(allMeasurements, modelOutputs);
            auto offsets = hierarchical.computeAnalyticalOffsets(allMeasurements, modelOutputs);
            hierarchical.applyOptimalScalings(scalings, modelOutputs);
            hierarchical.applyOptimalOffsets(offsets, modelOutputs);
            // TODO: what else needs to be scaled?
        }

        // TODO auto sigmas = hierarchical.computeAnalyticalOffsets(allMeasurements, modelOutputs);

        for(int dataIdx = 0; (unsigned) dataIdx < jobs.size(); ++dataIdx) {
            JobResultAmiciSimulation& result = jobResults[dataIdx];
            result.rdata->y = modelOutputs[dataIdx];
            std::vector<double> sigmas(allMeasurements[dataIdx].size(), NAN);
            result.rdata->llh = -parpe::computeNegLogLikelihood(allMeasurements[dataIdx], modelOutputs[dataIdx], sigmas);

            auto edata = dataProvider->getExperimentalDataForCondition(dataIdx);
            rw.saveSimulationResults(edata.get(), result.rdata.get(), dataIdx);
        }
        return 0;
    });


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
    auto model = dataProvider->getModel();
    auto solver = dataProvider->getSolver();
    JobAmiciSimulation<JobIdentifier> sim(solver.get(), model.get(), &path);
    sim.deserialize(buffer.data(), buffer.size());
    solver->setSensitivityOrder(amici::AMICI_SENSI_ORDER_NONE);

#if QUEUE_WORKER_H_VERBOSE >= 2
    int mpiRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    printf("[%d] Received work. ", mpiRank);
    fflush(stdout);
#endif

    // do work
    JobResultAmiciSimulation result = runSimulation(path, *solver, *model);

#if QUEUE_WORKER_H_VERBOSE >= 2
    printf("[%d] Work done. ", mpiRank);
    fflush(stdout);
#endif

    // pack & cleanup
    result.rdata->J = std::vector<double>();
    result.rdata->sigmay = std::vector<double>();
    result.rdata->ssigmay = std::vector<double>();
    result.rdata->sx0 = std::vector<double>();
    result.rdata->x = std::vector<double>();
    result.rdata->x0 = std::vector<double>();
    result.rdata->xdot = std::vector<double>();

    buffer = amici::serializeToStdVec<JobResultAmiciSimulation>(result);
}


JobResultAmiciSimulation StandaloneSimulator::runSimulation(JobIdentifier path,
                                                            amici::Solver& solver, amici::Model& model)
{
    // currently requires edata, since all condition specific parameters are set via edata
    auto edata = dataProvider->getExperimentalDataForCondition(path.idxConditions);

    auto rdata = amici::runAmiciSimulation(solver, edata.get(), model);

    RELEASE_ASSERT(rdata != nullptr, "");
    return JobResultAmiciSimulation(std::move(rdata), 0.0);
}



std::vector<double> getFinalParameters(std::string const& startIndex, H5::H5File const& file)
{
    auto lock = hdf5MutexGetLock();

    // find last iteration /multistarts/$/iteration/$/costFunParameters
    std::string iterationPath = std::string("/multistarts/") + startIndex + "/iteration/";
    int iteration = 0;
    while(hdf5GroupExists(file.getId(), (iterationPath + std::to_string(iteration)).c_str()) && hdf5DatasetExists(file.getId(), iterationPath + std::to_string(iteration) + "/costFunParameters")) {
        ++iteration;
    }
    --iteration; // last one did not exist

    auto bestPairLast = getFunctionEvaluationWithMinimalCost(
                iterationPath + std::to_string(iteration) + "/costFunCost",
                file);
    int costFunEvaluationIndex = bestPairLast.first;

    if(iteration > 0) {
        // If job got killed during line search, the final point of the previous iteration
        // might be better than any line search steps of the current iteration
        auto bestPairSecondLast = getFunctionEvaluationWithMinimalCost(
                    iterationPath + std::to_string(iteration - 1) + "/costFunCost",
                    file);
        if(bestPairSecondLast.second < bestPairLast.second) {
            --iteration;
            costFunEvaluationIndex = bestPairSecondLast.first;
        }
    }

    // get parameters of the selected function evaluation
    std::string parameterPath = iterationPath + std::to_string(iteration) + "/costFunParameters";
    H5::DataSet dataset = file.openDataSet(parameterPath);

    H5::DataSpace filespace = dataset.getSpace();
    int rank = filespace.getSimpleExtentNdims();
    RELEASE_ASSERT(rank == 2, "Rank mismatch");

    hsize_t dims[2];
    filespace.getSimpleExtentDims(dims, nullptr);
    int numParam = dims[0];

    std::vector<double> parameters(numParam);

    parpe::hdf5Read2DDoubleHyperslab(file.getId(), parameterPath.c_str(),
                                     numParam, 1, 0, costFunEvaluationIndex,
                                     parameters.data());

    /*
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
*/
    return parameters;
}

std::pair<int, double> getFunctionEvaluationWithMinimalCost(std::string const& datasetPath, H5::H5File const& file) {
    H5::DataSet dataset = file.openDataSet(datasetPath);

    H5::DataSpace filespace = dataset.getSpace();
    int rank = filespace.getSimpleExtentNdims();
    RELEASE_ASSERT(rank == 2, "Rank mismatch");

    hsize_t dims[2];
    filespace.getSimpleExtentDims(dims, nullptr);
    RELEASE_ASSERT(dims[0] == 1, "Dim1 should be 1");
    int numFunctionEvalations = dims[1];

    std::vector<double> cost(numFunctionEvalations, INFINITY);

    parpe::hdf5Read2DDoubleHyperslab(file.getId(), datasetPath.c_str(),
                                     1, numFunctionEvalations, 0, 0,
                                     cost.data());
    int minIndex = std::min_element(cost.begin(), cost.end()) - cost.begin();
    return std::pair<int, double>(minIndex, cost[minIndex]);
}

std::vector<std::vector<double>> getParameterTrajectory(std::string const& startIndex, H5::H5File const& file)
{
    auto lock = hdf5MutexGetLock();

    std::string parameterPath = std::string("/multistarts/") + startIndex + "/iterCostFunParameters";
    H5::DataSet dataset = file.openDataSet(parameterPath);

    H5::DataSpace filespace = dataset.getSpace();
    int rank = filespace.getSimpleExtentNdims();
    RELEASE_ASSERT(rank == 2, "Rank mismatch");

    hsize_t dims[2];
    filespace.getSimpleExtentDims(dims, nullptr);
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

int getNumStarts(H5::H5File const& file, std::string const& rootPath)  {
    auto o = parpe::OptimizationOptions::fromHDF5(file.getId(), rootPath + "/inputData/optimizationOptions");
    return o->numStarts;
}

int runFinalParameters(StandaloneSimulator &sim,
                       std::string const& inFileName,
                       std::string const& resultFileName,
                       std::string const& resultPath,
                       LoadBalancerMaster *loadBalancer) {

    H5::H5File file;
    {
        auto lock = hdf5MutexGetLock();
        file.openFile(inFileName, H5F_ACC_RDONLY);
    }
    int errors = 0;

    int numStarts = getNumStarts(file);
    for(int i = 0; i < numStarts; ++i) {
        std::cout<<"Running for start "<<i<<std::endl;
        try {
            auto parameters = parpe::getFinalParameters(std::to_string(i), file);
            std::string curResultPath = resultPath + "multistarts/" + std::to_string(i);
            errors += sim.run(resultFileName, curResultPath, parameters, loadBalancer, file);
        } catch (H5::FileIException e) {
            std::cerr<<"Exception during start " << i << " "<<e.getDetailMsg()<<std::endl;
            std::cerr<<"... skipping"<<std::endl;
        }
    }

    return errors;
}

int runAlongTrajectory(StandaloneSimulator &sim,
                       const std::string &inFileName,
                       const std::string &resultFileName,
                       const std::string &resultPath,
                       LoadBalancerMaster *loadBalancer)
{
    H5::H5File file;
    {
        auto lock = hdf5MutexGetLock();
        file.openFile(inFileName, H5F_ACC_RDONLY);
    }

    int errors = 0;

    for(int i = 0; i < getNumStarts(file); ++i) {
        try {

            auto parameters = getParameterTrajectory(std::to_string(i), file);

            for(int iter = 0; (unsigned) iter < parameters.size(); ++iter) {
                std::cout<<"Running for start "<<i<<" iter "<<iter<<std::endl;
                std::string curResultPath = resultPath + "multistarts/" + std::to_string(i) + "/iter/" + std::to_string(iter);

                errors += sim.run(resultFileName, curResultPath, parameters[iter], loadBalancer, file);
            }
        } catch (std::exception const& e) {
            std::cerr<<e.what()<<std::endl;
        }

    }

    return errors;
}

int runSimulationTasks(StandaloneSimulator &sim,
                       std::string const& simulationMode,
                       std::string const& inFileName,
                       std::string const& dataFilePath,
                       std::string const& resultFileName,
                       std::string const& resultPath,
                       LoadBalancerMaster *loadBalancer) {

    if(simulationMode == "--at-optimum") {
        return parpe::runFinalParameters(sim, inFileName, resultFileName, resultPath, loadBalancer);
    }

    if (simulationMode == "--along-trajectory") {
        return parpe::runAlongTrajectory(sim, inFileName, resultFileName, resultPath, loadBalancer);
    }

    return -1;
}

int runSimulator(MultiConditionDataProvider &dp,
                 const std::string &simulationMode,
                 const std::string &inFileName,
                 const std::string &dataFilePath,
                 const std::string &resultFileName,
                 const std::string &resultPath)
{
    parpe::StandaloneSimulator sim(&dp);
    int status = 0;
    int commSize = parpe::getMpiCommSize();

    if (commSize > 1) {
        if (parpe::getMpiRank() == 0) {
            parpe::LoadBalancerMaster loadBalancer;
            loadBalancer.run();
            status = runSimulationTasks(sim, simulationMode,
                                        inFileName, dataFilePath,
                                        resultFileName, resultPath, &loadBalancer);
            loadBalancer.terminate();
            loadBalancer.sendTerminationSignalToAllWorkers();
        } else {
            parpe::LoadBalancerWorker lbw;
            lbw.run([&sim](std::vector<char> &buffer, int jobId) {
                sim.messageHandler(buffer, jobId);
            });
        }
    } else {
        status = runSimulationTasks(sim, simulationMode,
                                    inFileName, dataFilePath,
                                    resultFileName, resultPath, nullptr);
    }

    return status;
}


} // namespace parpe
