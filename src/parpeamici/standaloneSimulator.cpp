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
                         std::vector<double> const& optimizationParameters,
                         LoadBalancerMaster* loadBalancer,
                         H5::H5File const& conditionFile,
                         std::string conditionFilePath)
{
    // std::cout<<"file: "<<resultFile<<" path: "<<resultPath<<"
    // lbm:"<<loadBalancer<<std::endl;

    SimulationResultWriter rw(resultFile, resultPath);
    rw.saveYMes = true;
    rw.saveYSim = true;
    rw.saveLlh = true;
    rw.save_parameters_ = true;
    rw.saveX = true;

    auto model = dataProvider->getModel();
    auto solver = dataProvider->getSolver();
    solver->setSensitivityOrder(amici::SensitivityOrder::none);

    std::vector<double> parameters = optimizationParameters;
    HierarchicalOptimizationWrapper hierarchical(nullptr, 0, 0);

    /* If the provided parameter vector is shorter than required, this means, we
     * got only the result of the outer problem of an hierarchical optimization
     * run, and thus, need to compute the inner optimal parameters here.
     */
    // TODO: check parameter names! until this is implemented, let's always
    // recompute output parameters
    bool needComputeAnalyticalParameters = true;
    // (parameters.size() !=
    // (unsigned)dataProvider->getNumOptimizationParameters());

    if (needComputeAnalyticalParameters) {
        if (hdf5GroupExists(conditionFile.getId(), "inputData"))
            // TODO: might not be the best place to have that here
            conditionFilePath += "/inputData";
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

        auto wrappedFun = std::make_unique<AmiciSummedGradientFunction>(
          dataProvider, loadBalancer, nullptr);
        wrappedFun->sendStates = true;

        hierarchical = HierarchicalOptimizationWrapper(
          std::move(wrappedFun),
          std::move(hierarchicalScalingReader),
          std::move(hierarchicalOffsetReader),
          std::move(hierarchicalSigmaReader),
          dataProvider->getNumberOfSimulationConditions(),
          model->nytrue,
          ErrorModel::normal);
        std::cout << "Need to compute analytical parameters: "
                  << conditionFilePath << "  "
                  << proportionalityFactorIndices.size()
                  << " parameters.size() == " << parameters.size()
                  << " ; hierarchical.numParameters() == "
                  << hierarchical.numParameters() << std::endl;
        Expects(parameters.size() == (unsigned)hierarchical.numParameters());

        // expand parameter vector
        auto scalingDummy = hierarchical.getDefaultScalingFactors();
        auto offsetDummy = hierarchical.getDefaultOffsetParameters();
        auto sigmaDummy = hierarchical.getDefaultSigmaParameters();
        parameters =
          spliceParameters(gsl::make_span(optimizationParameters.data(),
                                          optimizationParameters.size()),
                           proportionalityFactorIndices,
                           offsetParameterIndices,
                           sigmaParameterIndices,
                           scalingDummy,
                           offsetDummy,
                           sigmaDummy);
        // get outputs, scale
        // TODO need to pass aggregate function for writing
    } else {
        // is already the correct length
        // parameters = optimizationParameters;

        auto resultFileH5 = rw.reopenFile();
        hdf5EnsureGroupExists(resultFileH5.getId(), resultPath.c_str());
        auto lock = hdf5MutexGetLock();
        amici::hdf5::createAndWriteDouble1DDataset(
          resultFileH5, resultPath + "/problemParameters", parameters);
    }

    RELEASE_ASSERT(
        parameters.size() ==
            (unsigned)dataProvider->getNumOptimizationParameters(),
        "Size of supplied parameter vector does not match model dimensions.");

    rw.createDatasets(dataProvider->getNumberOfSimulationConditions());

    std::vector<int> dataIndices(
      dataProvider->getNumberOfSimulationConditions());
    std::iota(dataIndices.begin(), dataIndices.end(), 0);
    int errors = 0;
    std::cout << "Starting simulation. Number of conditions: "
              << dataProvider->getNumberOfSimulationConditions() << std::endl;
    auto jobFinished =
      [&](
        JobData* job,
        int /*dataIdx*/) { /* job finished */
                           // if we are running hierarchical optimization we
                           // need to wait until all jobs are finished
                           if (needComputeAnalyticalParameters)
                               return;

                           auto results = amici::deserializeFromChar<std::map<
                             int,
                             AmiciSimulationRunner::AmiciResultPackageSimple>>(
                             job->recvBuffer.data(), job->recvBuffer.size());
                           job->recvBuffer = std::vector<char>(); // free buffer

                           for (auto const& result : results) {
                               errors += result.second.status;
                               int conditionIdx = result.first;
                               auto edata =
                                 dataProvider->getExperimentalDataForCondition(
                                   conditionIdx);

                               rw.saveTimepoints(edata->getTimepoints(),
                                                 conditionIdx);
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
                               dataProvider->updateSimulationParametersAndScale(conditionIdx, parameters, *model);
                               rw.saveParameters(model->getParameters(),
                                                 conditionIdx);
                           }
      };

    auto allFinished = [&](std::vector<JobData>& jobs)
      -> int { /* all finished */
               if (!needComputeAnalyticalParameters)
                   return 0; // Work was already done in above function

               // must wait for all jobs to finish because of hierarchical
               // optimization and scaling factors
               std::vector<AmiciSimulationRunner::AmiciResultPackageSimple>
                 simulationResults(dataIndices.size());
               std::vector<std::vector<double>> modelOutputs(
                 dataIndices.size());
               std::vector<std::vector<double>> modelStates(
                   dataIndices.size());

               // collect all model outputs
               for (auto& job : jobs) {
                   auto results = amici::deserializeFromChar<
                     std::map<int,
                              AmiciSimulationRunner::AmiciResultPackageSimple>>(
                     job.recvBuffer.data(), job.recvBuffer.size());
                   job.recvBuffer = std::vector<char>(); // free buffer
                   for (auto& result : results) {
                       swap(simulationResults[result.first], result.second);
                       modelOutputs[result.first] =
                         simulationResults[result.first].modelOutput;
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
               auto fullSigmaMatrices = hierarchical.fun->getAllSigmas();
               if (!hierarchical.getSigmaParameterIndices().empty()) {
                   hierarchical.fillInAnalyticalSigmas(fullSigmaMatrices,
                                                       sigmas);
               }

               // save parameters
               parameters = spliceParameters(
                 parameters,
                 hierarchical.getProportionalityFactorIndices(),
                 hierarchical.getOffsetParameterIndices(),
                 hierarchical.getSigmaParameterIndices(),
                 scalings,
                 offsets,
                 sigmas);
               auto resultFileH5 = rw.reopenFile();
               hdf5EnsureGroupExists(resultFileH5.getId(), resultPath.c_str());
               {
                   auto lock = hdf5MutexGetLock();
                   amici::hdf5::createAndWriteDouble1DDataset(
                     resultFileH5, resultPath + "/problemParameters", parameters);
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
               // compute llh
               for (int conditionIdx = 0;
                    (unsigned)conditionIdx < simulationResults.size();
                    ++conditionIdx) {
                   double llh = -parpe::computeNegLogLikelihood(
                     allMeasurements[conditionIdx],
                     modelOutputs[conditionIdx],
                     fullSigmaMatrices[conditionIdx]);

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
                       conditionIdx, parameters, *model);
                   rw.saveParameters(model->getParameters(),
                                     conditionIdx);

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
        errors += simRunner.runDistributedMemory(loadBalancer,
                                                 maxSimulationsPerPackage);
    } else {
#endif
        errors +=
          simRunner.runSharedMemory([&](std::vector<char>& buffer, int jobId) {
              messageHandler(buffer, jobId);
          });
#ifdef PARPE_ENABLE_MPI
    }
#endif

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
            logger.logmessage(LOGLVL_ERROR, "[" + identifier + "] " + message);
        } else {
            logger.logmessage(LOGLVL_ERROR, message);
        }
    };
    amiciApp.warning = [&logger](std::string const& identifier,
                                 std::string const& message) {
        if (!identifier.empty()) {
            logger.logmessage(LOGLVL_WARNING,
                              "[" + identifier + "] " + message);
        } else {
            logger.logmessage(LOGLVL_WARNING, message);
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
        rdata->x,
        rdata->status
    };
}

std::vector<double>
getFinalParameters(std::string const& startIndex, H5::H5File const& file)
{
    auto lock = hdf5MutexGetLock();

    // find last iteration /multistarts/$/iteration/$/costFunParameters
    std::string iterationPath =
      std::string("/multistarts/") + startIndex + "/iteration/";
    int iteration = 0;
    while (
      hdf5GroupExists(file.getId(),
                      (iterationPath + std::to_string(iteration)).c_str()) &&
      hdf5DatasetExists(file,
                        iterationPath + std::to_string(iteration) +
                          "/costFunParameters")) {
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
      file.getId(), datasetPath.c_str(), 1, numFunctionEvalations, 0, 0, cost);
    int minIndex = std::min_element(cost.begin(), cost.end()) - cost.begin();
    return { minIndex, cost[minIndex] };
}

std::vector<std::vector<double>>
getParameterTrajectory(std::string const& startIndex, H5::H5File const& file)
{
    auto lock = hdf5MutexGetLock();

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
        parpe::hdf5Read2DDoubleHyperslab(file.getId(),
                                         parameterPath.c_str(),
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
      file.getId(), rootPath + "/inputData/optimizationOptions");
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
                   LoadBalancerMaster* loadBalancer)
{
    auto lock = hdf5MutexGetLock();
    H5::H5File parameterFile(parameterFileName, H5F_ACC_RDONLY);
    H5::H5File conditionFile(conditionFileName, H5F_ACC_RDONLY);
    lock.unlock();

    int errors = 0;

    int numStarts = getNumStarts(parameterFile);
    for (int i = 0; i < numStarts; ++i) {
        std::cout << "Running for start " << i << std::endl;
        try {
            auto parameters =
              parpe::getFinalParameters(std::to_string(i), parameterFile);
            auto outerParameters =
              getOuterParameters(parameters, parameterFile, parameterFilePath);

            std::string curResultPath =
              resultPath + "multistarts/" + std::to_string(i);


            errors += sim.run(resultFileName,
                              curResultPath,
                              outerParameters,
                              loadBalancer,
                              conditionFile,
                              conditionFilePath);
        } catch (H5::FileIException const& e) {
            std::cerr << "Exception during start " << i << " "
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
                   LoadBalancerMaster* loadBalancer)
{
    auto lock = hdf5MutexGetLock();
    H5::H5File parameterFile(parameterFileName, H5F_ACC_RDONLY);
    H5::H5File conditionFile(conditionFileName, H5F_ACC_RDONLY);
    lock.unlock();

    int errors = 0;

    for (int startIdx = 0; startIdx < getNumStarts(parameterFile); ++startIdx) {
        try {
            auto parameters =
              getParameterTrajectory(std::to_string(startIdx), parameterFile);

            for (int iter = 0; (unsigned)iter < parameters.size(); ++iter) {
                std::cout << "Running for start " << startIdx << " iter "
                          << iter << std::endl;
                std::string curResultPath = resultPath + "/multistarts/" +
                                            std::to_string(startIdx) +
                                            "/iter/" + std::to_string(iter);

                auto outerParameters = getOuterParameters(
                  parameters[iter], parameterFile, parameterFilePath);

                errors += sim.run(resultFileName,
                                  curResultPath,
                                  outerParameters,
                                  loadBalancer,
                                  conditionFile,
                                  conditionFilePath);
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
runSimulationTasks(StandaloneSimulator& sim,
                   std::string const& simulationMode,
                   std::string const& conditionFileName,
                   std::string const& conditionFilePath,
                   std::string const& parameterFileName,
                   std::string const& parameterFilePath,
                   std::string const& resultFileName,
                   std::string const& resultPath,
                   LoadBalancerMaster* loadBalancer)
{
    {
        // copy input data
        auto lock = hdf5MutexGetLock();
        H5::H5File conditionFile = hdf5OpenForReading(conditionFileName);
        H5::H5File resultFile = hdf5OpenForAppending(resultFileName);

        std::vector<std::string> datasetsToCopy {"/inputData"};
        for (auto const& datasetToCopy : datasetsToCopy) {
            auto source = conditionFilePath + datasetToCopy;
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
                                         loadBalancer);
    }

    if (simulationMode == "--along-trajectory") {
        return parpe::runAlongTrajectory(sim,
                                         conditionFileName,
                                         conditionFilePath,
                                         parameterFileName,
                                         parameterFilePath,
                                         resultFileName,
                                         resultPath,
                                         loadBalancer);
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
            status = runSimulationTasks(sim,
                                        simulationMode,
                                        conditionFileName,
                                        conditionFilePath,
                                        parameterFileName,
                                        parameterFilePath,
                                        resultFileName,
                                        resultPath,
                                        &loadBalancer);
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
                                    nullptr);
#ifdef PARPE_ENABLE_MPI
    }
#endif

    return status;
}

std::vector<double>
getOuterParameters(const std::vector<double>& fullParameters,
                   const H5::H5File& parameterFile,
                   const std::string& parameterPath)
{
    // auto options = OptimizationOptions::fromHDF5(parameterFile.getId(),
    // parameterPath + "/optimizationOptions");
    AnalyticalParameterHdf5Reader hierarchicalScalingReader(
      parameterFile,
      parameterPath + "/inputData/scalingParameterIndices",
      parameterPath + "/inputData/scalingParametersMapToObservables");
    AnalyticalParameterHdf5Reader hierarchicalOffsetReader(
      parameterFile,
      parameterPath + "/inputData/offsetParameterIndices",
      parameterPath + "/inputData/offsetParametersMapToObservables");
    AnalyticalParameterHdf5Reader hierarchicalSigmaReader(
      parameterFile,
      parameterPath + "/inputData/sigmaParameterIndices",
      parameterPath + "/inputData/sigmaParametersMapToObservables");

    auto proportionalityFactorIndices =
      hierarchicalScalingReader.getOptimizationParameterIndices();
    auto offsetParameterIndices =
      hierarchicalOffsetReader.getOptimizationParameterIndices();
    auto sigmaParameterIndices =
      hierarchicalSigmaReader.getOptimizationParameterIndices();

    auto combinedIndices = proportionalityFactorIndices;
    combinedIndices.insert(combinedIndices.end(),
                           offsetParameterIndices.begin(),
                           offsetParameterIndices.end());
    combinedIndices.insert(combinedIndices.end(),
                           sigmaParameterIndices.begin(),
                           sigmaParameterIndices.end());
    std::sort(combinedIndices.begin(), combinedIndices.end());

    std::vector<double> result(fullParameters.size() - combinedIndices.size());
    parpe::fillFilteredParams(fullParameters, combinedIndices, result);

    return result;
}

} // namespace parpe
