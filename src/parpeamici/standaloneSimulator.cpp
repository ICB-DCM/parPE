#include <parpeamici/standaloneSimulator.h>

#include <parpeamici/amiciSimulationRunner.h>
#include <parpeamici/simulationResultWriter.h>
#include <parpeloadbalancer/loadBalancerMaster.h>
#include <parpeoptimization/optimizationOptions.h>
#include <parpecommon/misc.h>
#include <parpeamici/hierarchicalOptimization.h>

#ifdef PARPE_ENABLE_MPI
#include <parpeloadbalancer/loadBalancerWorker.h>
#endif

#include <amici/hdf5.h>
#include <gsl/gsl-lite.hpp>

#include <iostream>

namespace parpe {

StandaloneSimulator::StandaloneSimulator(MultiConditionDataProvider *dp)
    : dataProvider(dp)
{
    if(auto env = std::getenv("PARPE_MAX_SIMULATIONS_PER_PACKAGE")) {
        maxSimulationsPerPackage = std::stoi(env);
    }
}

int StandaloneSimulator::run(const std::string& resultFile,
                             const std::string& resultPath,
                             std::vector<double> const& optimizationParameters,
                             LoadBalancerMaster *loadBalancer,
                             H5::H5File const& conditionFile,
                             std::string conditionFilePath)
{
    // std::cout<<"file: "<<resultFile<<" path: "<<resultPath<<" lbm:"<<loadBalancer<<std::endl;

    SimulationResultWriter rw(resultFile, resultPath);
    rw.saveYMes = true;
    rw.saveYSim = true;
    rw.saveLlh = true;

    auto model = dataProvider->getModel();
    auto solver = dataProvider->getSolver();
    solver->setSensitivityOrder(amici::SensitivityOrder::none);

    std::vector<double> parameters = optimizationParameters;
    HierarchicalOptimizationWrapper hierarchical(nullptr, 0, 0, 0);

    /* If the provided parameter vector is shorter than required, this means, we got only the result
     * of the outer problem of an hierarchical optimization run, and thus, need to compute the
     * inner optimal parameters here.
     */
    bool needComputeAnalyticalParameters = (parameters.size() != (unsigned)dataProvider->getNumOptimizationParameters());

    if(needComputeAnalyticalParameters) {
        if(hdf5GroupExists(conditionFile.getId(), "inputData"))
            // TODO: might not be the best place to have that here
            conditionFilePath += "/inputData";
        // TODO: get rid of that. we want fun.evaluate(), independently of hierarchical or not
        auto hierarchicalScalingReader = std::make_unique<AnalyticalParameterHdf5Reader>(
                    conditionFile,
                    conditionFilePath + "/scalingParameterIndices",
                    conditionFilePath + "/scalingParametersMapToObservables");
        auto hierarchicalOffsetReader = std::make_unique<AnalyticalParameterHdf5Reader>(
                    conditionFile,
                    conditionFilePath + "/offsetParameterIndices",
                    conditionFilePath + "/offsetParametersMapToObservables");
        auto hierarchicalSigmaReader = std::make_unique<AnalyticalParameterHdf5Reader>(
                    conditionFile,
                    conditionFilePath + "/sigmaParameterIndices",
                    conditionFilePath + "/sigmaParametersMapToObservables");
        auto proportionalityFactorIndices = hierarchicalScalingReader->getOptimizationParameterIndices();
        auto offsetParameterIndices = hierarchicalOffsetReader->getOptimizationParameterIndices();
        auto sigmaParameterIndices = hierarchicalSigmaReader->getOptimizationParameterIndices();

        auto wrappedFun = std::make_unique<AmiciSummedGradientFunction<int>>(dataProvider, loadBalancer, nullptr);

        hierarchical = HierarchicalOptimizationWrapper(std::move(wrappedFun),
                                                std::move(hierarchicalScalingReader),
                                                std::move(hierarchicalOffsetReader),
                                                std::move(hierarchicalSigmaReader),
                                                dataProvider->getNumberOfConditions(), model->nytrue, model->nt(),
                                                ErrorModel::normal);
        std::cout<<"Need to compute analytical parameters: "<<conditionFilePath<<"  "<<proportionalityFactorIndices.size()<<" parameters.size() == "<<parameters.size()<<" ; hierarchical.numParameters() == "<<hierarchical.numParameters()<<std::endl;
        RELEASE_ASSERT(parameters.size() == (unsigned) hierarchical.numParameters(), "");

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

    auto resultFileH5 = rw.reopenFile();
    hdf5EnsureGroupExists(resultFileH5.getId(), resultPath.c_str());
    amici::hdf5::createAndWriteDouble1DDataset(resultFileH5, resultPath + "/parameters", parameters.data(), parameters.size());

    RELEASE_ASSERT(parameters.size() == (unsigned)dataProvider->getNumOptimizationParameters(), "Size of supplied parameter vector does not match model dimensions.");

    rw.createDatasets(*model, dataProvider->getNumberOfConditions());

    std::vector<int> dataIndices(dataProvider->getNumberOfConditions());
    std::iota(dataIndices.begin(), dataIndices.end(), 0);
    int errors = 0;
    std::cout<<"Starting simulation. Number of conditions: " << dataProvider->getNumberOfConditions()<<std::endl;

    auto jobFinished = [&](JobData *job, int /*dataIdx*/) { /* job finished */
        // if we are running hierarchical optimization we need to wait until all jobs are finished
        if(needComputeAnalyticalParameters)
            return;

        auto results =
                amici::deserializeFromChar<
                std::map<int, AmiciSimulationRunner::AmiciResultPackageSimple>>(
                    job->recvBuffer.data(), job->recvBuffer.size());
        job->recvBuffer = std::vector<char>(); // free buffer

        for (auto const& result : results) {
            errors += result.second.status;
            int conditionIdx = result.first;
            auto edata = dataProvider->getExperimentalDataForCondition(conditionIdx);
            rw.saveMeasurements(edata->getObservedData(), edata->nt(), edata->nytrue(), conditionIdx);
            rw.saveModelOutputs(result.second.modelOutput,  model->nt(), model->nytrue, conditionIdx);
            rw.saveLikelihood(result.second.llh, conditionIdx);
        }
    };

    auto allFinished = [&](std::vector<JobData> &jobs) -> int { /* all finished */
        if(!needComputeAnalyticalParameters)
            return 0; // Work was already done in above function

        // must wait for all jobs to finish because of hierarchical optimization and scaling factors
        std::vector<AmiciSimulationRunner::AmiciResultPackageSimple> simulationResults(dataIndices.size());
        std::vector<std::vector<double> > modelOutputs(dataIndices.size());

        // collect all model outputs
        for(auto& job : jobs) {
            auto results = amici::deserializeFromChar<
                    std::map<int, AmiciSimulationRunner::AmiciResultPackageSimple>>(
                        job.recvBuffer.data(), job.recvBuffer.size());
            job.recvBuffer = std::vector<char>(); // free buffer
            for (auto& result : results) {
                swap(simulationResults[result.first], result.second);
                modelOutputs[result.first] = simulationResults[result.first].modelOutput;
            }
        }

        // TODO: redundant with hierarchicalOptimization.cpp
        //  compute scaling factors and offset parameters
        auto allMeasurements = dataProvider->getAllMeasurements();
        RELEASE_ASSERT(dataIndices.size() == allMeasurements.size(), "");

        auto scalings = hierarchical.computeAnalyticalScalings(allMeasurements, modelOutputs);
        auto offsets = hierarchical.computeAnalyticalOffsets(allMeasurements, modelOutputs);
        hierarchical.applyOptimalScalings(scalings, modelOutputs);
        hierarchical.applyOptimalOffsets(offsets, modelOutputs);
        auto sigmas = hierarchical.computeAnalyticalSigmas(allMeasurements, modelOutputs);
        auto fullSigmaMatrices = hierarchical.fun->getAllSigmas();
        if(!hierarchical.getSigmaParameterIndices().empty()) {
            hierarchical.fillInAnalyticalSigmas(fullSigmaMatrices, sigmas);
        }

        // compute llh
        for(int conditionIdx = 0; (unsigned) conditionIdx < simulationResults.size(); ++conditionIdx) {
            double llh = -parpe::computeNegLogLikelihood(
                        allMeasurements[conditionIdx],
                        modelOutputs[conditionIdx],
                        fullSigmaMatrices[conditionIdx]);
            auto edata = dataProvider->getExperimentalDataForCondition(conditionIdx);
            rw.saveMeasurements(edata->getObservedData(), edata->nt(), edata->nytrue(), conditionIdx);
            rw.saveModelOutputs(modelOutputs[conditionIdx],  model->nt(), model->nytrue, conditionIdx);
            rw.saveLikelihood(llh, conditionIdx);
        }
        return 0;
    };

    AmiciSimulationRunner simRunner(parameters,
                                    amici::SensitivityOrder::none,
                                    dataIndices,
                                    jobFinished,
                                    allFinished);

#ifdef PARPE_ENABLE_MPI
    if (loadBalancer && loadBalancer->isRunning()) {
        errors += simRunner.runDistributedMemory(loadBalancer, maxSimulationsPerPackage);
    } else {
#endif
        errors += simRunner.runSharedMemory(
                    [&](std::vector<char> &buffer, int jobId) {
                messageHandler(buffer, jobId);
    });
#ifdef PARPE_ENABLE_MPI
    }
#endif

    return errors;
}


void StandaloneSimulator::messageHandler(std::vector<char> &buffer, int  /*jobId*/)
{
    // TODO: pretty redundant with messageHandler in multiconditionproblem
    // unpack simulation job data
    auto model = dataProvider->getModel();
    auto solver = dataProvider->getSolver();
    auto sim = amici::deserializeFromChar<
            AmiciSimulationRunner::AmiciWorkPackageSimple>(buffer.data(), buffer.size());
    solver->setSensitivityOrder(sim.sensitivityOrder);

#if QUEUE_WORKER_H_VERBOSE >= 2
    int mpiRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    printf("[%d] Received work. ", mpiRank);
    fflush(stdout);
#endif

    std::map<int, AmiciSimulationRunner::AmiciResultPackageSimple> results;
    // run simulations for all condition indices
    for(auto conditionIndex: sim.conditionIndices) {
        dataProvider->updateSimulationParameters(conditionIndex, sim.optimizationParameters, *model);
        auto result = runSimulation(conditionIndex, *solver, *model);
        results[conditionIndex] = result;
    }

#if QUEUE_WORKER_H_VERBOSE >= 2
    printf("[%d] Work done. ", mpiRank);
    fflush(stdout);
#endif

    buffer = amici::serializeToStdVec(results);
}


AmiciSimulationRunner::AmiciResultPackageSimple StandaloneSimulator::runSimulation(int conditionIdx,
                                                            amici::Solver& solver, amici::Model& model)
{
    // currently requires edata, since all condition specific parameters are set via edata
    auto edata = dataProvider->getExperimentalDataForCondition(conditionIdx);

    auto rdata = amici::runAmiciSimulation(solver, edata.get(), model);

    RELEASE_ASSERT(rdata != nullptr, "");

    return AmiciSimulationRunner::AmiciResultPackageSimple {
        rdata->llh,
                NAN,
                (solver.getSensitivityOrder() > amici::SensitivityOrder::none) ? rdata->sllh : std::vector<double>(),
                rdata->y,
                rdata->status
    };
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
    return { minIndex, cost[minIndex] };
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
                       std::string const& conditionFileName,
                       std::string const& conditionFilePath,
                       std::string const& parameterFileName,
                       std::string const& parameterFilePath,
                       std::string const& resultFileName,
                       std::string const& resultPath,
                       LoadBalancerMaster *loadBalancer) {

    H5::H5File parameterFile;
    {
        auto lock = hdf5MutexGetLock();
        parameterFile.openFile(parameterFileName, H5F_ACC_RDONLY);
    }
    int errors = 0;

    int numStarts = getNumStarts(parameterFile);
    for(int i = 0; i < numStarts; ++i) {
        std::cout<<"Running for start "<<i<<std::endl;
        try {
            auto parameters = parpe::getFinalParameters(std::to_string(i), parameterFile);
            auto outerParameters = getOuterParameters(parameters, parameterFile, parameterFilePath);

            std::string curResultPath = resultPath + "multistarts/" + std::to_string(i);

            auto lock = hdf5MutexGetLock();
            H5::H5File conditionFile = hdf5OpenForReading(conditionFileName);
            lock.unlock();

            errors += sim.run(resultFileName, curResultPath, outerParameters, loadBalancer,
                              conditionFile, conditionFilePath);
        } catch (H5::FileIException const& e) {
            std::cerr<<"Exception during start " << i << " "<<e.getDetailMsg()<<std::endl;
            std::cerr<<"... skipping"<<std::endl;
        }
    }

    return errors;
}

int runAlongTrajectory(StandaloneSimulator &sim,
                       const std::string &conditionFileName,
                       const std::string &conditionFilePath,
                       const std::string &parameterFileName,
                       const std::string &parameterFilePath,
                       std::string const& resultFileName,
                       std::string const& resultPath,
                       LoadBalancerMaster *loadBalancer)
{
    H5::H5File parameterFile;
    {
        auto lock = hdf5MutexGetLock();
        parameterFile.openFile(parameterFileName, H5F_ACC_RDONLY);
    }

    int errors = 0;

    for(int i = 0; i < getNumStarts(parameterFile); ++i) {
        try {
            auto parameters = getParameterTrajectory(std::to_string(i), parameterFile);

            for(int iter = 0; (unsigned) iter < parameters.size(); ++iter) {
                std::cout<<"Running for start "<<i<<" iter "<<iter<<std::endl;
                std::string curResultPath = resultPath + "/multistarts/" + std::to_string(i) + "/iter/" + std::to_string(iter);

                auto lock = hdf5MutexGetLock();
                H5::H5File conditionFile = hdf5OpenForReading(conditionFileName);
                lock.unlock();

                auto outerParameters = getOuterParameters(parameters[iter], parameterFile, parameterFilePath);

                errors += sim.run(resultFileName, curResultPath,
                                  outerParameters, loadBalancer,
                                  conditionFile, conditionFilePath);
            }
        } catch (std::exception const& e) {
            std::cerr<<e.what()<<std::endl;
        }

    }

    return errors;
}

int runSimulationTasks(StandaloneSimulator &sim,
                       std::string const& simulationMode,
                       std::string const& conditionFileName,
                       std::string const& conditionFilePath,
                       std::string const& parameterFileName,
                       std::string const& parameterFilePath,
                       std::string const& resultFileName,
                       std::string const& resultPath,
                       LoadBalancerMaster *loadBalancer) {

    if(simulationMode == "--at-optimum") {
        return parpe::runFinalParameters(sim, conditionFileName, conditionFilePath,
                                         parameterFileName, parameterFilePath,
                                         resultFileName, resultPath, loadBalancer);
    }

    if (simulationMode == "--along-trajectory") {
        return parpe::runAlongTrajectory(sim, conditionFileName, conditionFilePath,
                                         parameterFileName, parameterFilePath,
                                         resultFileName, resultPath, loadBalancer);
    }

    return EXIT_FAILURE;
}

int runSimulator(MultiConditionDataProvider &dp,
                 std::string const& simulationMode,
                 std::string const& conditionFileName,
                 std::string const& conditionFilePath,
                 std::string const& parameterFileName,
                 std::string const& parameterFilePath,
                 std::string const& resultFileName,
                 std::string const& resultPath)
{
    parpe::StandaloneSimulator sim(&dp);
    int status = 0;

#ifdef PARPE_ENABLE_MPI
    int commSize = parpe::getMpiCommSize();
    if (commSize > 1) {
        if (parpe::getMpiRank() == 0) {
            parpe::LoadBalancerMaster loadBalancer;
            loadBalancer.run();
            status = runSimulationTasks(sim, simulationMode,
                                        conditionFileName, conditionFilePath,
                                        parameterFileName, parameterFilePath,
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
#endif
        status = runSimulationTasks(sim, simulationMode,
                                    conditionFileName, conditionFilePath,
                                    parameterFileName, parameterFilePath,
                                    resultFileName, resultPath, nullptr);
#ifdef PARPE_ENABLE_MPI
    }
#endif

    return status;
}

std::vector<double> getOuterParameters(const std::vector<double> &fullParameters,
                                       const H5::H5File &parameterFile,
                                       const std::string &parameterPath)
{
    //auto options = OptimizationOptions::fromHDF5(parameterFile.getId(), parameterPath + "/optimizationOptions");
    AnalyticalParameterHdf5Reader hierarchicalScalingReader(parameterFile,
                                                            parameterPath + "/inputData/scalingParameterIndices",
                                                            parameterPath + "/inputData/scalingParametersMapToObservables");
    AnalyticalParameterHdf5Reader hierarchicalOffsetReader(parameterFile,
                                                           parameterPath + "/inputData/offsetParameterIndices",
                                                           parameterPath + "/inputData/offsetParametersMapToObservables");
    AnalyticalParameterHdf5Reader hierarchicalSigmaReader(parameterFile,
                                                          parameterPath + "/inputData/sigmaParameterIndices",
                                                          parameterPath + "/inputData/sigmaParametersMapToObservables");

    auto proportionalityFactorIndices = hierarchicalScalingReader.getOptimizationParameterIndices();
    auto offsetParameterIndices = hierarchicalOffsetReader.getOptimizationParameterIndices();
    auto sigmaParameterIndices = hierarchicalSigmaReader.getOptimizationParameterIndices();

    auto combinedIndices = proportionalityFactorIndices;
    combinedIndices.insert(combinedIndices.end(), offsetParameterIndices.begin(), offsetParameterIndices.end());
    combinedIndices.insert(combinedIndices.end(), sigmaParameterIndices.begin(), sigmaParameterIndices.end());
    std::sort(combinedIndices.begin(), combinedIndices.end());

    std::vector<double> result(fullParameters.size() - combinedIndices.size());
    parpe::fillFilteredParams(fullParameters, combinedIndices, result);

    return result;
}


} // namespace parpe