#include <parpeamici/multiConditionProblem.h>

#include <parpecommon/logging.h>
#include <parpecommon/misc.h>

#include <parpeoptimization/optimizationOptions.h>
#include <parpeoptimization/optimizationResultWriter.h>

#include <parpeamici/hierarchicalOptimization.h>
#include <parpeamici/multiConditionDataProvider.h>
#include <parpeamici/amiciMisc.h>

#include <gsl/gsl-lite.hpp>

#include <cassert>
#include <cstring>
#include <ctime>
#include <numeric>
#include <utility>

#ifdef PARPE_ENABLE_MPI
#include <parpeloadbalancer/loadBalancerWorker.h>
#endif
#include <parpeloadbalancer/loadBalancerMaster.h>

#ifndef PARPE_ENABLE_MPI
// Workaround to allow building without MPI. Should be cleaned up.
using LoadBalancerMaster = int;
#endif

namespace parpe {

// For debugging:
// skip objective function evaluation completely
//#define NO_OBJ_FUN_EVAL

MultiConditionProblem::MultiConditionProblem(MultiConditionDataProvider *dp)
    : MultiConditionProblem(dp, nullptr, nullptr, nullptr) {}

MultiConditionProblem::MultiConditionProblem(
        MultiConditionDataProvider *dp,
        LoadBalancerMaster *loadBalancer,
        std::unique_ptr<Logger> logger,
        std::unique_ptr<OptimizationResultWriter> resultWriter)
    : dataProvider(dp),
      resultWriter(std::move(resultWriter))
{
    this->logger_ = std::move(logger);
    // run on all data
    std::vector<int> dataIndices(dataProvider->getNumberOfSimulationConditions());
    std::iota(dataIndices.begin(), dataIndices.end(), 0);

    cost_fun_ = std::make_unique<
            SummedGradientFunctionGradientFunctionAdapter<int>
            > (
                std::make_unique<AmiciSummedGradientFunction>(
                    dataProvider, loadBalancer, this->resultWriter.get()),
                dataIndices);

    if(auto hdp = dynamic_cast<MultiConditionDataProviderHDF5*>(dp)) {
        parametersMin.resize(dp->getNumOptimizationParameters());
        hdp->getOptimizationParametersLowerBounds(parametersMin);

        parametersMax.resize(dp->getNumOptimizationParameters());
        hdp->getOptimizationParametersUpperBounds(parametersMax);
    }

}

void MultiConditionProblem::fillParametersMin(gsl::span<double> buffer) const
{
    Expects(buffer.size() == parametersMin.size());
    std::copy(parametersMin.begin(), parametersMin.end(), buffer.begin());
}

void MultiConditionProblem::fillParametersMax(gsl::span<double> buffer) const
{
    Expects(buffer.size() == parametersMax.size());
    std::copy(parametersMax.begin(), parametersMax.end(), buffer.begin());
}

void MultiConditionProblem::fillInitialParameters(gsl::span<double> buffer) const
{
    if(!startingPoint.empty()) {
        Expects(buffer.size() == startingPoint.size());
        std::copy(startingPoint.begin(), startingPoint.end(), buffer.begin());
    } else {
        OptimizationProblem::fillInitialParameters(buffer);
    }
}


int MultiConditionProblem::earlyStopping() {
    bool stop = false;

    /* TODO evaluate objective function on test set and see if prediction
     * performance increases
     * costValidation <- validationProblem.evaluate();
     * costValidation.append()
     * if no decrease during last 3 rounds, return stop
     *
     * TODO: need to have current parameters; need to be supplied to
     * intermediate function;
     * need to from optimizer if in line search or actual step
     */

    // validationProblem->evaluateObjectiveFunction()

    return static_cast<int>(stop);
}


void MultiConditionProblem::setInitialParameters(
        std::vector<double> const& startingPoint)
{
    this->startingPoint = startingPoint;
}

void MultiConditionProblem::setParametersMin(
        std::vector<double> const& lowerBounds)
{
    parametersMin = lowerBounds;
}

void MultiConditionProblem::setParametersMax(
        std::vector<double> const& upperBounds)
{
    parametersMax = upperBounds;
}

std::unique_ptr<OptimizationReporter> MultiConditionProblem::getReporter() const
{

    return std::make_unique<OptimizationReporter>(
                cost_fun_.get(),
                std::make_unique<OptimizationResultWriter>(*resultWriter),
                std::make_unique<Logger>(*logger_));
}

std::vector<int> MultiConditionProblem::getTrainingData() const
{
    std::vector<int> dataIndices(
                dataProvider->getNumberOfSimulationConditions());
    std::iota(dataIndices.begin(), dataIndices.end(), 0);
    return dataIndices;
}

MultiConditionDataProvider *MultiConditionProblem::getDataProvider() {
    return dataProvider;
}

OptimizationResultWriter *MultiConditionProblem::getResultWriter() {
    return resultWriter.get();
}

MultiConditionProblemMultiStartOptimizationProblem
::MultiConditionProblemMultiStartOptimizationProblem(
        MultiConditionDataProviderHDF5 *dp,
        OptimizationOptions options,
        OptimizationResultWriter *resultWriter,
        LoadBalancerMaster *loadBalancer, std::unique_ptr<Logger> logger)
    : data_provider_(dp), options_(std::move(options)),
      result_writer_(resultWriter), load_balancer_(loadBalancer),
      logger_(std::move(logger))
{}

int MultiConditionProblemMultiStartOptimizationProblem::getNumberOfStarts() const { return options_.numStarts; }

bool MultiConditionProblemMultiStartOptimizationProblem::restartOnFailure() const { return options_.retryOptimization; }

std::unique_ptr<OptimizationProblem> MultiConditionProblemMultiStartOptimizationProblem::getLocalProblem(
        int multiStartIndex) const {
    // generate new OptimizationProblem with data from dp

    Expects(data_provider_ != nullptr);

    std::unique_ptr<MultiConditionProblem> problem;

    if (result_writer_) {
        problem = std::make_unique<MultiConditionProblem>(
                    data_provider_, load_balancer_,
                    logger_->getChild(std::string("o")
                                     + std::to_string(multiStartIndex)),
                    std::make_unique<OptimizationResultWriter>(*result_writer_));
        problem->getResultWriter()->setRootPath(
                    "/multistarts/" + std::to_string(multiStartIndex));
    } else {
        problem = std::make_unique<MultiConditionProblem>(
                    data_provider_, load_balancer_,
                    logger_->getChild(
                        std::string("o") + std::to_string(multiStartIndex)),
                    nullptr);
    }
    problem->setOptimizationOptions(options_);
    problem->setInitialParameters(
                parpe::OptimizationOptions::getStartingPoint(
                    data_provider_->getHdf5File(), multiStartIndex));

    if(options_.hierarchicalOptimization)
        return std::unique_ptr<OptimizationProblem>(
                    new parpe::HierarchicalOptimizationProblemWrapper(
                        std::move(problem), data_provider_));

    return problem;
}

void printSimulationResult(Logger *logger, int jobId,
                           amici::ReturnData const* rdata, double timeSeconds) {
    if(!rdata) {
        // This should not happen, but apparently we can't rely on AMICI always
        // returning some result object
        logger->logmessage(loglevel::error,
                           "AMICI simulation failed unexpectedly.");
        return;
    }

    bool with_sensi = rdata->sensi >= amici::SensitivityOrder::first;

    logger->logmessage(loglevel::debug, "Result for %d: %g (%d) (%d/%d/%.4fs%c)",
                       jobId, rdata->llh, rdata->status,
                       rdata->numsteps.empty()?-1:rdata->numsteps[rdata->numsteps.size() - 1],
                       rdata->numstepsB.empty()?-1:rdata->numstepsB[0],
                       timeSeconds,
                       with_sensi?'+':'-');


    // check for NaNs, only report first
    if (with_sensi) {
        for (int i = 0; i < rdata->np; ++i) {
            if (std::isnan(rdata->sllh[i])) {
                logger->logmessage(loglevel::debug,
                                   "Gradient contains NaN at %d", i);
                break;
            }

            if (std::isinf(rdata->sllh[i])) {
                logger->logmessage(loglevel::debug,
                                   "Gradient contains Inf at %d", i);
                break;
            }
        }
    }
}


void saveSimulation(const H5::H5File &file, std::string const& pathStr,
                    std::vector<double> const& parameters,
                    double llh, gsl::span<double const> gradient,
                    double timeElapsedInSeconds,
                    gsl::span<double const> /*states*/,
                    gsl::span<double const> /*stateSensi*/,
                    gsl::span<double const> /*outputs*/, int jobId,
                    int status, std::string const& label)
{
    // TODO replace by SimulationResultWriter
    const char *fullGroupPath = pathStr.c_str();

    [[maybe_unused]] auto lock = hdf5MutexGetLock();

    hdf5CreateOrExtendAndWriteToDouble2DArray(
        file, fullGroupPath, "simulationLogLikelihood",
                gsl::make_span<double>(&llh, 1));

    hdf5CreateOrExtendAndWriteToInt2DArray(
                file, fullGroupPath, "jobId",
                gsl::make_span<const int>(&jobId, 1));

    if (!gradient.empty()) {
        hdf5CreateOrExtendAndWriteToDouble2DArray(
            file, fullGroupPath, "simulationLogLikelihoodGradient",
                    gradient);
    } else if(!parameters.empty()) {
        double dummyGradient[parameters.size()];
        std::fill_n(dummyGradient, parameters.size(), NAN);
        hdf5CreateOrExtendAndWriteToDouble2DArray(
            file, fullGroupPath, "simulationLogLikelihoodGradient",
            gsl::make_span<const double>(dummyGradient, parameters.size()));
    }

    if (!parameters.empty())
        hdf5CreateOrExtendAndWriteToDouble2DArray(
            file, fullGroupPath, "simulationParameters", parameters);

    hdf5CreateOrExtendAndWriteToDouble2DArray(
                file, fullGroupPath, "simulationWallTimeInSec",
                gsl::make_span<const double>(&timeElapsedInSeconds, 1));

    // TODO: This was broken by allowing different numbers of timepoints
    // for different simulation conditions. Vector lengths now may differ and
    // lead to crashes.
    //    if (!states.empty())
    //        hdf5CreateOrExtendAndWriteToDouble2DArray(
    //            file.getId(), fullGroupPath, "simulationStates", states);

    //    if (!outputs.empty())
    //        hdf5CreateOrExtendAndWriteToDouble2DArray(
    //            file.getId(), fullGroupPath, "simulationObservables", outputs);

    //    if (!stateSensi.empty())
    //        hdf5CreateOrExtendAndWriteToDouble3DArray(
    //            file.getId(), fullGroupPath, "simulationStateSensitivities", stateSensi,
    //            stateSensi.size() / parameters.size(), parameters.size());

    hdf5CreateOrExtendAndWriteToInt2DArray(
                file, fullGroupPath, "simulationStatus",
                gsl::make_span<const int>(&status, 1));

    hdf5CreateOrExtendAndWriteToString1DArray(file, fullGroupPath,
                                           "simulationLabel", label);

    file.flush(H5F_SCOPE_LOCAL);

}

AmiciSimulationRunner::AmiciResultPackageSimple runAndLogSimulation(
        amici::Solver const &solverTemplate,
        amici::Model &model,
        int conditionIdx,
        int jobId,
        MultiConditionDataProvider const *dataProvider,
        OptimizationResultWriter *resultWriter,
        bool logLineSearch,
        Logger* logger,
        bool sendStates)
{
    // wall time  on worker for current simulation
    WallTimer simulationTimer;

    /* Get ExpData with measurement and fixed parameters. Other model parameters
     * and sensitivity options have been set already */
    auto edata = dataProvider->getExperimentalDataForCondition(conditionIdx);

    // TODO: extract class to handle tolerance relaxation

    /*
     * In case of simulation failure, try rerunning with a
     * `errorRelaxation`-fold higher error tolerance for a total of
     * `maxNumTrials` times (including the initial attempt).
     */
    constexpr int defaultMaxNumTrials = 6;
    constexpr double defaultErrorRelaxation = 1e2;

    int maxNumTrials = defaultMaxNumTrials;
    double errorRelaxation = defaultErrorRelaxation;

    // Set via environment variables?
    if(auto env = std::getenv("PARPE_NUM_SIMULATION_TRIALS")) {
        maxNumTrials = std::stoi(env);
    }
    if(auto env =
            std::getenv("PARPE_INTEGRATION_TOLERANCE_RELAXATION_FACTOR")) {
        errorRelaxation = std::stod(env);
    }

    std::unique_ptr<amici::ReturnData> rdata;

    // redirect AMICI output to parPE logging
    amici::AmiciApplication amiciApp;
    amiciApp.error = [logger](
            std::string const& identifier,
            std::string const& message){
        if(!identifier.empty()) {
            logger->logmessage(loglevel::error, "[" + identifier + "] " + message);
        } else {
            logger->logmessage(loglevel::error, message);
        }
    };
    amiciApp.warning = [logger](
            std::string const& identifier,
            std::string const& message){
        if(!identifier.empty()) {
            logger->logmessage(loglevel::warning,
                               "[" + identifier + "] " + message);
        } else {
            logger->logmessage(loglevel::warning, message);
        }
    };
    model.app = &amiciApp; // TODO: may dangle need to unset on exit

    for(int trial = 1; trial <= maxNumTrials; ++trial) {
        /* It is currently not safe to reuse solver if an exception has
         * occurred,so clone every time */
        auto solver = std::unique_ptr<amici::Solver>(solverTemplate.clone());
        solver->app = &amiciApp;
        if (!sendStates) {
            /* If we don't need the states, we can save memory here.
             * For current optimizers we only need the likelihood. For
             * hierarchical optimization we need the model outputs. Here, we
             * don't know about this, but for know it seems safe to use
             * amici::RDataReporting::likelihood if sensitivities are requested
             * and RDataReporting::residuals otherwise
             */

            if(solver->getSensitivityOrder() >= amici::SensitivityOrder::first
                && solver->getSensitivityMethod()
                       == amici::SensitivityMethod::adjoint) {
                solver->setReturnDataReportingMode(amici::RDataReporting::likelihood);
            } else {
                // unset sensitivity method, because `residuals` is not allowed
                // with `adjoint`, independent of sensitivity order
                solver->setSensitivityMethod(amici::SensitivityMethod::none);
                solver->setReturnDataReportingMode(amici::RDataReporting::residuals);
            }

        }

        if(trial - 1 == maxNumTrials) {
            logger->logmessage(loglevel::error,
                               "Simulation trial %d/%d failed. Giving up.",
                               trial, maxNumTrials);
            break;
        }

        if(rdata) {
            /* something went wrong in the previous simulation. until we have
             * better exception handling, we check those fields to deduce where
             * the error occurred
             */
            bool forwardFailed = std::isnan(rdata->llh);
            bool backwardFailed = solver->getSensitivityOrder() >= amici::SensitivityOrder::first
                    && solver->getSensitivityMethod() == amici::SensitivityMethod::adjoint
                    && !rdata->sllh.empty() && std::isnan(rdata->sllh[0]);

            // relax respective tolerances
            if(forwardFailed) {
                solver->setAbsoluteTolerance(
                            std::pow(errorRelaxation, trial - 1)
                            * solver->getAbsoluteTolerance());
                solver->setRelativeTolerance(
                            std::pow(errorRelaxation, trial - 1)
                            * solver->getRelativeTolerance());
            } else if (backwardFailed) {
                solver->setAbsoluteToleranceQuadratures(
                            std::pow(errorRelaxation, trial - 1)
                            * solver->getAbsoluteToleranceQuadratures());
                solver->setRelativeToleranceQuadratures(
                            std::pow(errorRelaxation, trial - 1)
                            * solver->getRelativeToleranceQuadratures());
                solver->setAbsoluteToleranceB(
                            std::pow(errorRelaxation, trial - 1)
                            * solver->getAbsoluteToleranceB());
                solver->setRelativeToleranceB(
                            std::pow(errorRelaxation, trial - 1)
                            * solver->getRelativeToleranceB());
            }

            logger->logmessage(
                        loglevel::warning,
                        "Error during simulation (try %d/%d), "
                        "retrying with relaxed error tolerances (*= %g): "
                        "abs: %g rel: %g quadAbs: %g quadRel: %g "
                        "abs_asa: %g, rel_asa: %g",
                        trial - 1, maxNumTrials, errorRelaxation,
                        solver->getAbsoluteTolerance(),
                        solver->getRelativeTolerance(),
                        solver->getAbsoluteToleranceQuadratures(),
                        solver->getRelativeToleranceQuadratures(),
                        solver->getAbsoluteToleranceB(),
                        solver->getRelativeToleranceB());
        }

        try {
            rdata = amiciApp.runAmiciSimulation(*solver, edata.get(), model);
        } catch (std::exception const& e) {
            std::cerr<<e.what()<<std::endl;
            std::string status = "-";
            if(rdata) {
                status = std::to_string(rdata->status);
            }
            logger->logmessage(
                        loglevel::warning, "Error during simulation: %s (%s)",
                e.what(), status.c_str());
        }

        if(rdata && rdata->status == amici::AMICI_SUCCESS)
            break;
    }
    double timeSeconds = simulationTimer.getTotal();

    printSimulationResult(logger, jobId, rdata.get(), timeSeconds);

    if (resultWriter && rdata && (solverTemplate.getSensitivityOrder()
                         > amici::SensitivityOrder::none || logLineSearch)) {
        saveSimulation(resultWriter->getH5File(), resultWriter->getRootPath(),
                       model.getParameters(), rdata->llh, rdata->sllh,
                       timeSeconds, rdata->x, rdata->sx, rdata->y,
                       jobId, rdata->status, logger->getPrefix());
    }
    if(rdata) {
        return AmiciSimulationRunner::AmiciResultPackageSimple {
            rdata->llh,
            timeSeconds,
            (solverTemplate.getSensitivityOrder()
             > amici::SensitivityOrder::none)
                ? rdata->sllh : std::vector<double>(),
            rdata->y, rdata->sigmay,
            sendStates ? rdata->x : std::vector<double>(),
            rdata->status
        };
    }

    // AMICI failed expectedly and did not return anything
    return AmiciSimulationRunner::AmiciResultPackageSimple {
        NAN,
        timeSeconds,
        (solverTemplate.getSensitivityOrder()
         > amici::SensitivityOrder::none)
            ? std::vector<double>(model.nplist(), NAN) : std::vector<double>(),
        std::vector<double>(model.nytrue, NAN),
        std::vector<double>(model.nytrue, NAN),
        sendStates ? std::vector<double>(model.nx_rdata, NAN) : std::vector<double>(),
        amici::AMICI_UNRECOVERABLE_ERROR
    };

}

FunctionEvaluationStatus getModelOutputsAndSigmas(
        MultiConditionDataProvider *dataProvider,
        LoadBalancerMaster *loadBalancer,
        int maxSimulationsPerPackage,
        OptimizationResultWriter *resultWriter,
        bool logLineSearch,
        gsl::span<const double> parameters,
        std::vector<std::vector<double> > &modelOutputs,
        std::vector<std::vector<double> > &modelSigmas,
        Logger *logger, double * /*cpuTime*/, bool sendStates)
{
    int errors = 0;

    std::vector<int> dataIndices(dataProvider->getNumberOfSimulationConditions());
    std::iota(dataIndices.begin(), dataIndices.end(), 0);

    modelOutputs.resize(dataIndices.size());
    auto parameterVector = std::vector<double>(parameters.begin(),
                                               parameters.end());
    auto jobFinished = [&errors, &modelOutputs, &modelSigmas](JobData *job, int /*dataIdx*/) {
        // deserialize
        auto results =
                amici::deserializeFromChar<AmiciSummedGradientFunction::ResultMap> (
                    job->recvBuffer.data(), job->recvBuffer.size());
        std::vector<char>().swap(job->recvBuffer); // free buffer

        for (auto const& result : results) {
            errors += result.second.status;
            modelOutputs[result.first] = result.second.modelOutput;
            modelSigmas[result.first] = result.second.modelSigmas;
        }
    };

    AmiciSimulationRunner simRunner(parameterVector,
                                    amici::SensitivityOrder::none,
                                    dataIndices,
                                    jobFinished,
                                    nullptr /* aggregate */,
                                    logger?logger->getPrefix():"");


#ifdef PARPE_ENABLE_MPI
    if (loadBalancer && loadBalancer->isRunning()) {
        errors += simRunner.runDistributedMemory(loadBalancer,
                                                 maxSimulationsPerPackage);
    } else {
#endif
        errors += simRunner.runSharedMemory(
                    [&dataProvider, &resultWriter, logLineSearch, &sendStates](std::vector<char> &buffer, int jobId) {
                messageHandler(dataProvider, resultWriter, logLineSearch,
                               buffer, jobId, sendStates);
    });
#ifdef PARPE_ENABLE_MPI
    }
#endif
    return errors == 0 ? functionEvaluationSuccess
                       : functionEvaluationFailure;
}



void messageHandler(MultiConditionDataProvider *dataProvider,
                    OptimizationResultWriter *resultWriter,
                    bool logLineSearch,
                    std::vector<char> &buffer, int jobId,
                    bool sendStates) {

#if QUEUE_WORKER_H_VERBOSE >= 2
    int mpiRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    printf("[%d] Received work. ", mpiRank);
    fflush(stdout);
#endif

    auto solver = dataProvider->getSolver();
    auto model = dataProvider->getModel();

    // unpack simulation job data
    auto workPackage = amici::deserializeFromChar<AmiciSummedGradientFunction::WorkPackage>(
                buffer.data(), buffer.size());

    solver->setSensitivityOrder(workPackage.sensitivityOrder);

    AmiciSummedGradientFunction::ResultMap results;
    // run simulations for all condition indices
    for(auto conditionIdx: workPackage.conditionIndices) {
        dataProvider->updateSimulationParametersAndScale(
                    conditionIdx,
                    workPackage.optimizationParameters,
                    *model);
        Logger logger(workPackage.logPrefix
                      + "c" + std::to_string(conditionIdx));
        auto result = runAndLogSimulation(
                    *solver, *model, conditionIdx, jobId, dataProvider,
                    resultWriter, logLineSearch, &logger, sendStates);
        results[conditionIdx] = result;
    }

#if QUEUE_WORKER_H_VERBOSE >= 2
    printf("[%d] Work done. ", mpiRank);
    fflush(stdout);
#endif
    // serialize to output buffer
    buffer = amici::serializeToStdVec(results);
}

AmiciSummedGradientFunction::AmiciSummedGradientFunction(
        MultiConditionDataProvider *dataProvider,
        LoadBalancerMaster *loadBalancer,
        OptimizationResultWriter *resultWriter)
    : dataProvider(dataProvider),
      loadBalancer(loadBalancer),
      model(dataProvider->getModel()),
      solver(dataProvider->getSolver()),
      solverOriginal(solver->clone()),
      resultWriter(resultWriter)
{
    if(auto env = std::getenv("PARPE_LOG_SIMULATIONS")) {
        logLineSearch = env[0] == '1';
    }

    if(auto env = std::getenv("PARPE_MAX_SIMULATIONS_PER_PACKAGE")) {
        maxSimulationsPerPackage = std::stoi(env);
    }

    if(auto env =
            std::getenv("PARPE_MAX_GRADIENT_SIMULATIONS_PER_PACKAGE")) {
        maxGradientSimulationsPerPackage = std::stoi(env);
    }
}

FunctionEvaluationStatus AmiciSummedGradientFunction::evaluate(
        gsl::span<const double> parameters, int dataset,
        double &fval, gsl::span<double> gradient, Logger *logger,
        double *cpuTime) const
{
    std::vector<int> datasets(1);
    datasets.at(0) = dataset;
    return evaluate(parameters, datasets, fval, gradient, logger, cpuTime);
}

FunctionEvaluationStatus AmiciSummedGradientFunction::evaluate(
        gsl::span<const double> parameters, std::vector<int> datasets,
        double &fval, gsl::span<double> gradient, Logger *logger,
        double *cpuTime) const
{
#ifdef NO_OBJ_FUN_EVAL
    if (objectiveFunctionGradient)
        std::fill(objectiveFunctionGradient,
                  objectiveFunctionGradient + numOptimizationParameters_,
                  0);
    *objectiveFunctionValue = 1;
    return 0;
#endif

    setSensitivityOptions(!gradient.empty());
    fval = 0.0;
    if (!gradient.empty())
        std::fill(gradient.begin(), gradient.end(), 0.0);

    int errors = runSimulations(parameters, fval, gradient, datasets,
                                logger, cpuTime);

    if (errors || !std::isfinite(fval)) {
        fval = std::numeric_limits<double>::infinity();
        return functionEvaluationFailure;
    }

    return functionEvaluationSuccess;
}

int AmiciSummedGradientFunction::numParameters() const
{
    return dataProvider->getNumOptimizationParameters();
}

std::vector<std::string> AmiciSummedGradientFunction::getParameterIds() const
{
    return dataProvider->getProblemParameterIds();
}

FunctionEvaluationStatus AmiciSummedGradientFunction::getModelOutputsAndSigmas(
    gsl::span<const double> parameters,
    std::vector<std::vector<double> > &modelOutputs,
    std::vector<std::vector<double> > &modelSigmas,
    Logger *logger, double *cpuTime) const
{
    return parpe::getModelOutputsAndSigmas(
        dataProvider, loadBalancer, maxSimulationsPerPackage, resultWriter,
        logLineSearch, parameters, modelOutputs, modelSigmas,
        logger, cpuTime, sendStates);
}


std::vector<std::vector<double> > AmiciSummedGradientFunction::getAllMeasurements() const {
    return dataProvider->getAllMeasurements();
}

void AmiciSummedGradientFunction::messageHandler(std::vector<char> &buffer, int jobId) const {
    parpe::messageHandler(dataProvider, resultWriter, logLineSearch, buffer,
                          jobId, sendStates);
}

amici::ParameterScaling AmiciSummedGradientFunction::getParameterScaling(
        int parameterIndex) const
{
    // parameterIndex is optimization parameter index,
    // not necessarily model parameter index!
    return dataProvider->getParameterScaleOpt(parameterIndex);
}

int AmiciSummedGradientFunction::runSimulations(
        gsl::span<const double> optimizationParameters,
        double &nllh, gsl::span<double> objectiveFunctionGradient,
        std::vector<int> const& dataIndices, Logger *logger, double *cpuTime) const {

    int errors = 0;

    auto parameterVector = std::vector<double>(
                optimizationParameters.begin(),
                optimizationParameters.end());
    double simulationTimeSec = 0.0;

    AmiciSimulationRunner simRunner(
                parameterVector,
                !objectiveFunctionGradient.empty()
                ? amici::SensitivityOrder::first
                : amici::SensitivityOrder::none,
                dataIndices,
                [&nllh, &objectiveFunctionGradient, &simulationTimeSec,
         &optimizationParameters, &errors, this](JobData *job, int /*jobIdx*/) {
            errors += this->aggregateLikelihood(*job,
                                      nllh,
                                      objectiveFunctionGradient,
                                      simulationTimeSec,
                                      optimizationParameters);
    }, nullptr,  logger?logger->getPrefix():"");

#ifdef PARPE_ENABLE_MPI
    if (loadBalancer && loadBalancer->isRunning()) {
        // When running simulations (without gradient),
        // send more simulations to each worker
        // to reduce communication overhead
        errors += simRunner.runDistributedMemory(
                    loadBalancer,
                    !objectiveFunctionGradient.empty()
                    ? maxGradientSimulationsPerPackage
                    : maxSimulationsPerPackage);
    } else {
#endif
        errors += simRunner.runSharedMemory(
                    [this](std::vector<char> &buffer, int jobId) {
                this->messageHandler(buffer, jobId);
    });
#ifdef PARPE_ENABLE_MPI
    }
#endif
    if(cpuTime)
        *cpuTime = simulationTimeSec;

    return errors;
}

int AmiciSummedGradientFunction::aggregateLikelihood(
        JobData &data, double &negLogLikelihood,
        gsl::span<double> negLogLikelihoodGradient,
        double &simulationTimeInS,
        gsl::span<const double> optimizationParameters) const
{
    int errors = 0;

    // deserialize
    auto results =
            amici::deserializeFromChar<ResultMap> (
                data.recvBuffer.data(), data.recvBuffer.size());
    std::vector<char>().swap(data.recvBuffer); // free buffer


    for (auto const& result : results) {
        int conditionIdx;
        ResultPackage resultPackage;
        std::tie(conditionIdx, resultPackage) = result;

        errors += resultPackage.status != amici::AMICI_SUCCESS;

        // sum up
        negLogLikelihood -= resultPackage.llh;
        simulationTimeInS += resultPackage.simulationTimeSeconds;

        if (!negLogLikelihoodGradient.empty()) {
            std::vector<double> p(model->np());
            auto scaleSim = dataProvider->getParameterScaleSim(conditionIdx);
            auto scaleOpt = dataProvider->getParameterScaleOpt();

            dataProvider->mapAndSetOptimizationToSimulationVariables(
                        conditionIdx, optimizationParameters, p, scaleOpt,
                        scaleSim);
            addSimulationGradientToObjectiveFunctionGradient(
                        conditionIdx, resultPackage.gradient,
                        negLogLikelihoodGradient, p);
        }
    }
    return errors;
}

void AmiciSummedGradientFunction::addSimulationGradientToObjectiveFunctionGradient(
        int conditionIdx, gsl::span<const double> simulationGradient,
        gsl::span<double> objectiveFunctionGradient,
        gsl::span<const double> simulationParameters) const {
    dataProvider->mapSimulationToOptimizationGradientAddMultiply(
                conditionIdx, simulationGradient,
                objectiveFunctionGradient, simulationParameters, -1.0);
}

void AmiciSummedGradientFunction::setSensitivityOptions(bool sensiRequired) const {
    // sensitivities requested?
    if (sensiRequired) {
        solver->setSensitivityOrder(solverOriginal->getSensitivityOrder());
        solver->setSensitivityMethod(solverOriginal->getSensitivityMethod());
    } else {
        solver->setSensitivityOrder(amici::SensitivityOrder::none);
        solver->setSensitivityMethod(amici::SensitivityMethod::none);
    }
}


} // namespace parpe
