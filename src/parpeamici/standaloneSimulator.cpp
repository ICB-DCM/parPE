#include <parpeamici/standaloneSimulator.h>

#include <parpeamici/amiciSimulationRunner.h>
#include <parpeamici/hierarchicalOptimization.h>
#include <parpeamici/multiConditionDataProvider.h>
#include <parpeamici/simulationResultWriter.h>
#include <parpecommon/misc.h>
#include <parpeloadbalancer/loadBalancerMaster.h>
#include <parpeoptimization/optimizationOptions.h>

#ifdef PARPE_ENABLE_MPI
#include <parpeloadbalancer/loadBalancerWorker.h>
#endif

#include <amici/hdf5.h>
#include <gsl/gsl-lite.hpp>

#include <iostream>

namespace parpe {

StandaloneSimulator::StandaloneSimulator(MultiConditionDataProvider* dp)
    : dataProvider(dp)
{
    if (auto env = std::getenv("PARPE_MAX_SIMULATIONS_PER_PACKAGE")) {
        maxSimulationsPerPackage = std::stoi(env);
    }
}

int
StandaloneSimulator::run(const std::string& resultFile,
                         const std::string& resultPath,
                         std::map<std::string, double>& optimizationParameters,
                         LoadBalancerMaster* loadBalancer,
                         H5::H5File const& conditionFile,
                         std::string conditionFilePath,
                         bool computeInnerParameters)
{
    // Data to save
    SimulationResultWriter rw(resultFile, resultPath);
    rw.saveYMes = true;
    rw.saveYSim = true;
    rw.saveLlh = true;
    rw.save_parameters_ = true;
    rw.saveX = true;

    // Setup model and solver
    auto model = dataProvider->getModel();
    auto solver = dataProvider->getSolver();
    solver->setSensitivityOrder(amici::SensitivityOrder::none);

    auto wrappedFun = std::make_unique<AmiciSummedGradientFunction>(
        dataProvider, loadBalancer, nullptr);
    wrappedFun->sendStates = true;

    AmiciSimulationRunner::callbackJobFinishedType jobFinished;
    AmiciSimulationRunner::callbackAllFinishedType allFinished;

    std::vector<int> dataIndices(
        dataProvider->getNumberOfSimulationConditions());
    std::iota(dataIndices.begin(), dataIndices.end(), 0);

    std::vector<double> parameterValues;
    int errors = 0;

    rw.createDatasets(dataProvider->getNumberOfSimulationConditions());

    // Write IDs
    auto resultFileH5 = rw.reopenFile();
    hdf5EnsureGroupExists(resultFileH5.getId(), resultPath.c_str());
    {
        [[maybe_unused]] auto lock = hdf5MutexGetLock();
        hdf5Write1dStringDataset(resultFileH5,
                                 resultPath,
                                 "stateIds",
                                 model->getStateIds());
        hdf5Write1dStringDataset(resultFileH5,
                                 resultPath,
                                 "observableIds",
                                 model->getObservableIds());
        hdf5Write1dStringDataset(resultFileH5,
                                 resultPath,
                                 "parameterIds",
                                 model->getParameterIds());
    }

    std::cout << "Starting simulation. Number of conditions: "
              << dataProvider->getNumberOfSimulationConditions() << std::endl;
    std::unique_ptr<AmiciSimulationRunner> simRunner;

    if(computeInnerParameters) {
        // hierarchical optimization
        HierarchicalOptimizationWrapper hierarchical(nullptr, 0, 0);

        // TODO: get rid of that. we want fun.evaluate(), independently of
        // hierarchical or not
        auto hierarchicalScalingReader =
            std::make_unique<AnalyticalParameterHdf5Reader>(
                conditionFile,
                conditionFilePath + "/scalingParameterIndices",
                conditionFilePath + "/scalingParametersMapToObservables");
        auto hierarchicalOffsetReader =
            std::make_unique<AnalyticalParameterHdf5Reader>(
                conditionFile,
                conditionFilePath + "/offsetParameterIndices",
                conditionFilePath + "/offsetParametersMapToObservables");
        auto hierarchicalSigmaReader =
            std::make_unique<AnalyticalParameterHdf5Reader>(
                conditionFile,
                conditionFilePath + "/sigmaParameterIndices",
                conditionFilePath + "/sigmaParametersMapToObservables");
        auto proportionalityFactorIndices =
            hierarchicalScalingReader->getOptimizationParameterIndices();
        auto offsetParameterIndices =
            hierarchicalOffsetReader->getOptimizationParameterIndices();
        auto sigmaParameterIndices =
            hierarchicalSigmaReader->getOptimizationParameterIndices();

        hierarchical = HierarchicalOptimizationWrapper(
            wrappedFun.get(),
            std::move(hierarchicalScalingReader),
            std::move(hierarchicalOffsetReader),
            std::move(hierarchicalSigmaReader),
            dataProvider->getNumberOfSimulationConditions(),
            model->nytrue,
            ErrorModel::normal);

        // Collect parameter values
        auto outerParameterNames = hierarchical.getParameterIds();
        std::vector<double> outerParameters(outerParameterNames.size());
        for(int i = 0; i < static_cast<int>(outerParameterNames.size()); ++i) {
            outerParameters[i] = optimizationParameters[outerParameterNames[i]];
        }
        Expects(outerParameters.size() == (unsigned)hierarchical.numParameters());

        // expand parameter vector
        auto scalingDummy = hierarchical.getDefaultScalingFactors();
        auto offsetDummy = hierarchical.getDefaultOffsetParameters();
        auto sigmaDummy = hierarchical.getDefaultSigmaParameters();
        parameterValues = spliceParameters(
            outerParameters,
            proportionalityFactorIndices,
            offsetParameterIndices,
            sigmaParameterIndices,
            scalingDummy,
            offsetDummy,
            sigmaDummy);

        allFinished = [this, &dataIndices, &hierarchical, &parameterValues,
                       &outerParameters, &resultFileH5, &resultPath, &model,
                       &rw](std::vector<JobData>& jobs)
        {
            /* all finished */
            // must wait for all jobs to finish because of hierarchical
            // optimization and scaling factors
            std::vector<AmiciSimulationRunner::AmiciResultPackageSimple>
                simulationResults(dataIndices.size());
            std::vector<std::vector<double>> modelOutputs(
                dataIndices.size());
            std::vector<std::vector<double>> modelSigmas(
                dataIndices.size());
            std::vector<std::vector<double>> modelStates(
                dataIndices.size());

            // collect all model outputs
            for (auto& job : jobs) {
                auto results = amici::deserializeFromChar<
                    std::map<int,
                             AmiciSimulationRunner::AmiciResultPackageSimple>>(
                    job.recvBuffer.data(), job.recvBuffer.size());
                std::vector<char>().swap(job.recvBuffer); // free buffer
                for (auto& result : results) {
                    swap(simulationResults[result.first], result.second);
                    modelOutputs[result.first] =
                        simulationResults[result.first].modelOutput;
                    modelSigmas[result.first] =
                        simulationResults[result.first].modelSigmas;
                    modelStates[result.first] =
                        simulationResults[result.first].modelStates;

                }
            }

            // TODO: redundant with hierarchicalOptimization.cpp
            //  compute scaling factors and offset parameters
            auto allMeasurements = dataProvider->getAllMeasurements();
            Expects(dataIndices.size() == allMeasurements.size());

            auto scalings = hierarchical.computeAnalyticalScalings(
                allMeasurements, modelOutputs);
            auto offsets = hierarchical.computeAnalyticalOffsets(
                allMeasurements, modelOutputs);
            hierarchical.applyOptimalScalings(scalings, modelOutputs);
            hierarchical.applyOptimalOffsets(offsets, modelOutputs);
            auto sigmas = hierarchical.computeAnalyticalSigmas(
                allMeasurements, modelOutputs);
            if (!hierarchical.getSigmaParameterIndices().empty()) {
                hierarchical.fillInAnalyticalSigmas(modelSigmas,
                                                    sigmas);
            }

            // save parameters
            parameterValues = spliceParameters(
                outerParameters,
                hierarchical.getProportionalityFactorIndices(),
                hierarchical.getOffsetParameterIndices(),
                hierarchical.getSigmaParameterIndices(),
                scalings,
                offsets,
                sigmas);
            {
                [[maybe_unused]] auto lock = hdf5MutexGetLock();
                amici::hdf5::createAndWriteDouble1DDataset(
                    resultFileH5, resultPath + "/problemParameters", parameterValues);
            }

            // compute llh
            for (int conditionIdx = 0;
                 (unsigned)conditionIdx < simulationResults.size();
                 ++conditionIdx) {
                double llh = -parpe::computeNegLogLikelihood(
                    allMeasurements[conditionIdx],
                    modelOutputs[conditionIdx],
                    modelSigmas[conditionIdx]);

                auto edata = dataProvider->getExperimentalDataForCondition(
                    conditionIdx);
                rw.saveTimepoints(edata->getTimepoints(), conditionIdx);
                if(!modelStates[conditionIdx].empty()) {
                    rw.saveStates(modelStates[conditionIdx], edata->nt(),
                                  model->nx_rdata, conditionIdx);
                }
                rw.saveMeasurements(edata->getObservedData(),
                                    edata->nt(),
                                    edata->nytrue(),
                                    conditionIdx);
                rw.saveModelOutputs(modelOutputs[conditionIdx],
                                    edata->nt(),
                                    model->nytrue,
                                    conditionIdx);
                rw.saveLikelihood(llh, conditionIdx);

                // to save simulation parameters
                dataProvider->updateSimulationParametersAndScale(
                    conditionIdx, parameterValues, *model);
                rw.saveParameters(model->getParameters(),
                                  conditionIdx);

            }
            return 0;
        };

        simRunner = std::make_unique<AmiciSimulationRunner>(
            parameterValues,
            amici::SensitivityOrder::none,
            dataIndices,
            jobFinished,
            allFinished);

#ifdef PARPE_ENABLE_MPI
        if (loadBalancer && loadBalancer->isRunning()) {
            errors += simRunner->runDistributedMemory(
                loadBalancer, maxSimulationsPerPackage);
        } else {
#endif
            errors +=
                simRunner->runSharedMemory([this](std::vector<char>& buffer, int jobId) {
                    messageHandler(buffer, jobId);
                });
#ifdef PARPE_ENABLE_MPI
        }
#endif
    } else {
        // Collect parameter values
        auto parameterNames = wrappedFun->getParameterIds();
        parameterValues.resize(parameterNames.size());
        for(int i = 0; i < static_cast<int>(parameterNames.size()); ++i) {
            parameterValues[i] = optimizationParameters[parameterNames[i]];
        }
        Expects(parameterValues.size() == (unsigned)wrappedFun->numParameters());

        {
            [[maybe_unused]] auto lock = hdf5MutexGetLock();
            amici::hdf5::createAndWriteDouble1DDataset(
                resultFileH5, resultPath + "/problemParameters", parameterValues);
        }
        jobFinished =
            [&](
                JobData* job,
                int /*dataIdx*/) {
                /* job finished */
                auto results = amici::deserializeFromChar<std::map<
                    int,
                    AmiciSimulationRunner::AmiciResultPackageSimple>>(
                    job->recvBuffer.data(), job->recvBuffer.size());
                std::vector<char>().swap(job->recvBuffer); // free buffer

                for (auto const& result : results) {
                    errors += result.second.status;
                    int conditionIdx = result.first;
                    auto edata =
                        dataProvider->getExperimentalDataForCondition(
                            conditionIdx);

                    rw.saveTimepoints(edata->getTimepoints(),
                                      conditionIdx);
                    if(!result.second.modelStates.empty()) {
                        rw.saveStates(result.second.modelStates, edata->nt(),
                                      model->nx_rdata, conditionIdx);
                    }
                    rw.saveMeasurements(edata->getObservedData(),
                                        edata->nt(),
                                        edata->nytrue(),
                                        conditionIdx);
                    rw.saveModelOutputs(result.second.modelOutput,
                                        edata->nt(),
                                        model->nytrue,
                                        conditionIdx);
                    rw.saveLikelihood(result.second.llh,
                                      conditionIdx);

                    // to save simulation parameters
                    dataProvider->updateSimulationParametersAndScale(
                        conditionIdx, parameterValues, *model);
                    rw.saveParameters(model->getParameters(), conditionIdx);
                }
            };

        simRunner = std::make_unique<AmiciSimulationRunner>(
            parameterValues,
            amici::SensitivityOrder::none,
            dataIndices,
            jobFinished,
            allFinished);

#ifdef PARPE_ENABLE_MPI
        if (loadBalancer && loadBalancer->isRunning()) {
            errors += simRunner->runDistributedMemory(
                loadBalancer, maxSimulationsPerPackage);
        } else {
#endif
            errors +=
                simRunner->runSharedMemory([this](std::vector<char>& buffer, int jobId) {
                    messageHandler(buffer, jobId);
                });
#ifdef PARPE_ENABLE_MPI
        }
#endif
    }

    return errors;
}

void
StandaloneSimulator::messageHandler(std::vector<char>& buffer, int /*jobId*/)
{
    // TODO: pretty redundant with messageHandler in multiconditionproblem
    // unpack simulation job data
    auto model = dataProvider->getModel();
    auto solver = dataProvider->getSolver();
    auto sim =
        amici::deserializeFromChar<AmiciSimulationRunner::AmiciWorkPackageSimple>(
            buffer.data(), buffer.size());
    solver->setSensitivityOrder(sim.sensitivityOrder);

#if QUEUE_WORKER_H_VERBOSE >= 2
    int mpiRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    printf("[%d] Received work. ", mpiRank);
    fflush(stdout);
#endif

    std::map<int, AmiciSimulationRunner::AmiciResultPackageSimple> results;
    // run simulations for all condition indices
    for (auto conditionIndex : sim.conditionIndices) {
        dataProvider->updateSimulationParametersAndScale(
            conditionIndex, sim.optimizationParameters, *model);
        auto result = runSimulation(conditionIndex, *solver, *model);
        results[conditionIndex] = result;
    }

#if QUEUE_WORKER_H_VERBOSE >= 2
    printf("[%d] Work done. ", mpiRank);
    fflush(stdout);
#endif

    buffer = amici::serializeToStdVec(results);
}

AmiciSimulationRunner::AmiciResultPackageSimple
StandaloneSimulator::runSimulation(int conditionIdx,
                                   amici::Solver& solver,
                                   amici::Model& model)
{
    // currently requires edata, since all condition specific parameters are set
    // via edata
    auto edata = dataProvider->getExperimentalDataForCondition(conditionIdx);

    // redirect AMICI output to parPE logging
    Logger logger("c" + std::to_string(conditionIdx));
    amici::AmiciApplication amiciApp;
    amiciApp.error = [&logger](std::string const& identifier,
                               std::string const& message) {
        if (!identifier.empty()) {
            logger.logmessage(loglevel::error, "[" + identifier + "] " + message);
        } else {
            logger.logmessage(loglevel::error, message);
        }
    };
    amiciApp.warning = [&logger](std::string const& identifier,
                                 std::string const& message) {
        if (!identifier.empty()) {
            logger.logmessage(loglevel::warning,
                              "[" + identifier + "] " + message);
        } else {
            logger.logmessage(loglevel::warning, message);
        }
    };
    model.app = &amiciApp;  // TODO: may dangle need to unset on exit
    solver.app = &amiciApp; // TODO: may dangle need to unset on exit

    auto rdata = amiciApp.runAmiciSimulation(solver, edata.get(), model);

    Expects(rdata != nullptr);

    return AmiciSimulationRunner::AmiciResultPackageSimple{
        rdata->llh,
        NAN,
        (solver.getSensitivityOrder() > amici::SensitivityOrder::none)
            ? rdata->sllh
            : std::vector<double>(),
        rdata->y,
        rdata->sigmay,
        rdata->x,
        rdata->status
    };
}

std::vector<double>
getFinalParameters(std::string const& startIndex, H5::H5File const& file)
{
    [[maybe_unused]] auto lock = hdf5MutexGetLock();

    // find last iteration /multistarts/$/iteration/$/costFunParameters
    std::string iterationPath =
        std::string("/multistarts/") + startIndex + "/iteration/";
    int iteration = 0;
    while (
        hdf5GroupExists(file, iterationPath + std::to_string(iteration)) &&
        file.nameExists(iterationPath + std::to_string(iteration)
                        + "/costFunParameters")) {
        ++iteration;
    }
    --iteration; // last one did not exist

    auto bestPairLast = getFunctionEvaluationWithMinimalCost(
        iterationPath + std::to_string(iteration) + "/costFunCost", file);
    int costFunEvaluationIndex = bestPairLast.first;

    if (iteration > 0) {
        // If job got killed during line search, the final point of the previous
        // iteration might be better than any line search steps of the current
        // iteration
        auto bestPairSecondLast = getFunctionEvaluationWithMinimalCost(
            iterationPath + std::to_string(iteration - 1) + "/costFunCost", file);
        if (bestPairSecondLast.second < bestPairLast.second) {
            --iteration;
            costFunEvaluationIndex = bestPairSecondLast.first;
        }
    }

    // get parameters of the selected function evaluation
    std::string parameterPath =
        iterationPath + std::to_string(iteration) + "/costFunParameters";
    H5::DataSet dataset = file.openDataSet(parameterPath);

    H5::DataSpace filespace = dataset.getSpace();
    int rank = filespace.getSimpleExtentNdims();
    RELEASE_ASSERT(rank == 2, "Rank mismatch");

    hsize_t dims[2];
    filespace.getSimpleExtentDims(dims, nullptr);
    int numParam = dims[0];

    std::vector<double> parameters(numParam);

    parpe::hdf5Read2DDoubleHyperslab(file.getId(),
                                     parameterPath.c_str(),
                                     numParam,
                                     1,
                                     0,
                                     costFunEvaluationIndex,
                                     parameters);

    /*
    // read from last iteration (last column in
    /multistarts/$/iterCostFunParameters) std::string parameterPath =
    std::string("/multistarts/") + startIndex + "/iterCostFunParameters";
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

std::pair<int, double>
getFunctionEvaluationWithMinimalCost(std::string const& datasetPath,
                                     H5::H5File const& file)
{
    [[maybe_unused]] auto lock = hdf5MutexGetLock();

    H5::DataSet dataset = file.openDataSet(datasetPath);

    H5::DataSpace filespace = dataset.getSpace();
    int rank = filespace.getSimpleExtentNdims();
    RELEASE_ASSERT(rank == 2, "Rank mismatch");

    hsize_t dims[2];
    filespace.getSimpleExtentDims(dims, nullptr);
    RELEASE_ASSERT(dims[0] == 1, "Dim1 should be 1");
    int numFunctionEvalations = dims[1];

    std::vector<double> cost(numFunctionEvalations, INFINITY);

    parpe::hdf5Read2DDoubleHyperslab(
        file, datasetPath, 1, numFunctionEvalations, 0, 0, cost);
    int minIndex = std::min_element(cost.begin(), cost.end()) - cost.begin();
    return { minIndex, cost[minIndex] };
}

std::vector<std::vector<double>>
getParameterTrajectory(std::string const& startIndex, H5::H5File const& file)
{
    [[maybe_unused]] auto lock = hdf5MutexGetLock();

    std::string parameterPath =
        std::string("/multistarts/") + startIndex + "/iterCostFunParameters";
    H5::DataSet dataset = file.openDataSet(parameterPath);

    H5::DataSpace filespace = dataset.getSpace();
    int rank = filespace.getSimpleExtentNdims();
    RELEASE_ASSERT(rank == 2, "Rank mismatch");

    hsize_t dims[2];
    filespace.getSimpleExtentDims(dims, nullptr);
    int numIter = dims[1];
    int numParam = dims[0];

    std::vector<std::vector<double>> parameters(numIter);

    for (int iter = 0; iter < numIter; ++iter) {
        parameters[iter] = std::vector<double>(numParam);
        parpe::hdf5Read2DDoubleHyperslab(file,
                                         parameterPath,
                                         numParam,
                                         1,
                                         0,
                                         iter,
                                         parameters[iter]);
    }

    return parameters;
}

int
getNumStarts(H5::H5File const& file, std::string const& rootPath)
{
    auto o = parpe::OptimizationOptions::fromHDF5(
        file, rootPath + "/optimizationOptions");
    return o->numStarts;
}

int
runFinalParameters(StandaloneSimulator& sim,
                   std::string const& conditionFileName,
                   std::string const& conditionFilePath,
                   std::string const& parameterFileName,
                   std::string const& parameterFilePath,
                   std::string const& resultFileName,
                   std::string const& resultPath,
                   LoadBalancerMaster* loadBalancer, bool computeInnerParameters)
{
    [[maybe_unused]] auto lock = hdf5MutexGetLock();
    H5::H5File parameterFile(parameterFileName, H5F_ACC_RDONLY);
    H5::H5File conditionFile(conditionFileName, H5F_ACC_RDONLY);
    std::vector<std::string> parameterNames;
    if(hdf5GroupExists(parameterFile,
                        parameterFilePath + "/parameters/parameterNames")){
        parameterNames = hdf5Read1dStringDataset(
            parameterFile, parameterFilePath + "/parameters/parameterNames");
    } else {
        parameterNames = hdf5Read1dStringDataset(
            parameterFile, parameterFilePath + "/inputData/parameters/parameterNames");
    }
    lock.unlock();

    int errors = 0;

    int numStarts = getNumStarts(parameterFile);
    for (int iStart = 0; iStart < numStarts; ++iStart) {
        std::cout << "Running for start " << iStart << std::endl;
        try {
            auto parameterValues =
                parpe::getFinalParameters(std::to_string(iStart), parameterFile);
            Expects(parameterValues.size() == parameterNames.size());
            std::map<std::string, double> parameters;
            for(int i = 0; i < static_cast<int>(parameters.size()); ++i)
                parameters[parameterNames[i]] = parameterValues[i];

            std::string curResultPath =
                resultPath + "multistarts/" + std::to_string(iStart);

            errors += sim.run(resultFileName,
                              curResultPath,
                              parameters,
                              loadBalancer,
                              conditionFile,
                              conditionFilePath,
                              computeInnerParameters);
        } catch (H5::FileIException const& e) {
            std::cerr << "Exception during start " << iStart << " "
                      << e.getDetailMsg() << std::endl;
            std::cerr << "... skipping" << std::endl;
        }
    }

    // lock for destruction of H5Files
    // FIXME: won't lock if an unhandled exception occurs
    lock.lock();

    return errors;
}

int
runAlongTrajectory(StandaloneSimulator& sim,
                   const std::string& conditionFileName,
                   const std::string& conditionFilePath,
                   const std::string& parameterFileName,
                   const std::string& parameterFilePath,
                   std::string const& resultFileName,
                   std::string const& resultPath,
                   LoadBalancerMaster* loadBalancer, bool computeInnerParameters)
{
    [[maybe_unused]] auto lock = hdf5MutexGetLock();
    H5::H5File parameterFile(parameterFileName, H5F_ACC_RDONLY);
    H5::H5File conditionFile(conditionFileName, H5F_ACC_RDONLY);
    std::vector<std::string> parameterNames;
    if(hdf5GroupExists(parameterFile,
                        parameterFilePath + "/parameters/parameterNames")){
        parameterNames = hdf5Read1dStringDataset(
            parameterFile, parameterFilePath + "/parameters/parameterNames");
    } else {
        parameterNames = hdf5Read1dStringDataset(
            parameterFile, parameterFilePath + "/inputData/parameters/parameterNames");
    }

    lock.unlock();

    int errors = 0;

    for (int startIdx = 0; startIdx < getNumStarts(parameterFile); ++startIdx) {
        try {
            auto parameterTrajectory =
                getParameterTrajectory(std::to_string(startIdx), parameterFile);

            for (int iter = 0; (unsigned)iter < parameterTrajectory.size(); ++iter) {
                std::cout << "Running for start " << startIdx << " iter "
                          << iter << std::endl;
                std::string curResultPath = resultPath + "/multistarts/" +
                                            std::to_string(startIdx) +
                                            "/iter/" + std::to_string(iter);
                auto const& parameterValues = parameterTrajectory[iter];

                Expects(parameterValues.size() == parameterNames.size());
                std::map<std::string, double> parameters;
                for(int i = 0; i < static_cast<int>(parameters.size()); ++i)
                    parameters[parameterNames[i]] = parameterValues[i];

                errors += sim.run(resultFileName,
                                  curResultPath,
                                  parameters,
                                  loadBalancer,
                                  conditionFile,
                                  conditionFilePath,
                                  computeInnerParameters);
            }
        } catch (std::exception const& e) {
            std::cerr << e.what() << std::endl;
        }
    }

    // lock for destruction of H5Files
    // FIXME: won't lock if an unhandled exception occurs
    lock.lock();

    return errors;
}


int
runNominalParameters(StandaloneSimulator& sim,
                     std::string const& conditionFileName,
                     std::string const& conditionFilePath,
                     std::string const& parameterFileName,
                     std::string const& parameterFilePath,
                     std::string const& resultFileName,
                     std::string const& resultPath,
                     LoadBalancerMaster* loadBalancer,
                     bool computeInnerParameters)
{
    [[maybe_unused]] auto lock = hdf5MutexGetLock();
    H5::H5File parameterFile(parameterFileName, H5F_ACC_RDONLY);
    H5::H5File conditionFile(conditionFileName, H5F_ACC_RDONLY);

    int errors = 0;

    std::cout << "Running for nominal parameter from "
              <<parameterFileName<<" "
              <<parameterFilePath + "/parameters/nominalValues"<< std::endl;
    // Read nominal parameters
    auto parameterValues = amici::hdf5::getDoubleDataset1D(
        parameterFile, parameterFilePath + "/parameters/nominalValues");
    auto parameterNames = hdf5Read1dStringDataset(
        parameterFile, parameterFilePath + "/parameters/parameterNames");
    lock.unlock();

    Expects(parameterValues.size() == parameterNames.size());
    std::map<std::string, double> parameters;
    for(int i = 0; i < static_cast<int>(parameterValues.size()); ++i)
        parameters[parameterNames[i]] = parameterValues[i];

    std::string curResultPath = resultPath + "nominal/";


    errors += sim.run(resultFileName,
                      curResultPath,
                      parameters,
                      loadBalancer,
                      conditionFile,
                      conditionFilePath,
                      computeInnerParameters);

    // lock for destruction of H5Files
    // FIXME: won't lock if an unhandled exception occurs
    lock.lock();

    return errors;
}

int
runSimulationTasks(StandaloneSimulator& sim,
                   std::string const& simulationMode,
                   std::string const& conditionFileName,
                   std::string const& conditionFilePath,
                   std::string const& parameterFileName,
                   std::string const& parameterFilePath,
                   std::string const& resultFileName,
                   std::string const& resultPath,
                   LoadBalancerMaster* loadBalancer,
                   bool computeInnerParameters)
{
    {
        std::cout<<"Running "<<simulationMode<<" for \n\tconditions from "
                  <<conditionFileName<<":"<<conditionFilePath
                  <<" and \n\tparameters from "
                  <<parameterFileName<<":"<<parameterFilePath<<"\n\t> "
                  <<resultFileName<<":"<<resultPath<<std::endl;
        // copy input data
        [[maybe_unused]] auto lock = hdf5MutexGetLock();
        H5::H5File conditionFile = hdf5OpenForReading(conditionFileName);
        H5::H5File resultFile = hdf5OpenForAppending(resultFileName);

        // TODO: this may not always be present. decide elsewhere what to copy
        std::vector<std::string> datasetsToCopy {"/inputData"};
        for (auto const& datasetToCopy : datasetsToCopy) {
            auto source = conditionFilePath + datasetToCopy;

            if(!conditionFile.exists(source)) {
                continue;
            }

            auto dest = resultPath + "/" + datasetToCopy;
            H5Ocopy(conditionFile.getId(), source.c_str(),
                    resultFile.getId(), dest.c_str(),
                    H5P_DEFAULT, H5P_DEFAULT);
        }
    }

    if (simulationMode == "--at-optimum") {
        return parpe::runFinalParameters(sim,
                                         conditionFileName,
                                         conditionFilePath,
                                         parameterFileName,
                                         parameterFilePath,
                                         resultFileName,
                                         resultPath,
                                         loadBalancer,
                                         computeInnerParameters);
    }

    if (simulationMode == "--along-trajectory") {
        return parpe::runAlongTrajectory(sim,
                                         conditionFileName,
                                         conditionFilePath,
                                         parameterFileName,
                                         parameterFilePath,
                                         resultFileName,
                                         resultPath,
                                         loadBalancer,
                                         computeInnerParameters);
    }

    if (simulationMode == "--nominal") {
        return parpe::runNominalParameters(sim,
                                           conditionFileName,
                                           conditionFilePath,
                                           parameterFileName,
                                           parameterFilePath,
                                           resultFileName,
                                           resultPath,
                                           loadBalancer,
                                           computeInnerParameters);

    }

    return EXIT_FAILURE;
}

int
runSimulator(MultiConditionDataProvider& dp,
             std::string const& simulationMode,
             std::string const& conditionFileName,
             std::string const& conditionFilePath,
             std::string const& parameterFileName,
             std::string const& parameterFilePath,
             std::string const& resultFileName,
             std::string const& resultPath,
             bool computeInnerParameters)
{
    parpe::StandaloneSimulator sim(&dp);
    int status = 0;

#ifdef PARPE_ENABLE_MPI
    int commSize = parpe::getMpiCommSize();
    if (commSize > 1) {
        if (parpe::getMpiRank() == 0) {
            parpe::LoadBalancerMaster loadBalancer;
            loadBalancer.run();
            status = runSimulationTasks(sim,
                                        simulationMode,
                                        conditionFileName,
                                        conditionFilePath,
                                        parameterFileName,
                                        parameterFilePath,
                                        resultFileName,
                                        resultPath,
                                        &loadBalancer,
                                        computeInnerParameters);
            loadBalancer.terminate();
            loadBalancer.sendTerminationSignalToAllWorkers();
        } else {
            parpe::LoadBalancerWorker lbw;
            lbw.run([&sim](std::vector<char>& buffer, int jobId) {
                sim.messageHandler(buffer, jobId);
            });
        }
    } else {
#endif
        status = runSimulationTasks(sim,
                                    simulationMode,
                                    conditionFileName,
                                    conditionFilePath,
                                    parameterFileName,
                                    parameterFilePath,
                                    resultFileName,
                                    resultPath,
                                    nullptr,
                                    computeInnerParameters);
#ifdef PARPE_ENABLE_MPI
    }
#endif

    return status;
}


} // namespace parpe
