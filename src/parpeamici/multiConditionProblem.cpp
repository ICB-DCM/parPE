#include <parpeamici/multiConditionProblem.h>

#include <parpeamici/steadystateSimulator.h>
#include <parpeoptimization/optimizationOptions.h>
#include <parpecommon/logging.h>
#include <parpecommon/misc.h>
#include <parpeamici/hierarchicalOptimization.h>
#include <parpeamici/multiConditionDataProvider.h>
#include <parpeoptimization/optimizationResultWriter.h>

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
    this->logger = std::move(logger);
    // run on all data
    std::vector<int> dataIndices(dataProvider->getNumberOfConditions());
    std::iota(dataIndices.begin(), dataIndices.end(), 0);

    costFun = std::make_unique<
            SummedGradientFunctionGradientFunctionAdapter<int>
            > (
                std::make_unique<AmiciSummedGradientFunction>(
                    dataProvider, loadBalancer, this->resultWriter.get()),
                dataIndices);

    if(auto hdp = dynamic_cast<MultiConditionDataProviderHDF5*>(dp)) {
        parametersMin.resize(dp->getNumOptimizationParameters());
        hdp->getOptimizationParametersLowerBounds(parametersMin.data());

        parametersMax.resize(dp->getNumOptimizationParameters());
        hdp->getOptimizationParametersUpperBounds(parametersMax.data());
    }

}

void MultiConditionProblem::fillParametersMin(gsl::span<double> buffer) const
{
    RELEASE_ASSERT(buffer.size() == parametersMin.size(), "");
    std::copy(parametersMin.begin(), parametersMin.end(), buffer.begin());
}

void MultiConditionProblem::fillParametersMax(gsl::span<double> buffer) const
{
    RELEASE_ASSERT(buffer.size() == parametersMax.size(), "");
    std::copy(parametersMax.begin(), parametersMax.end(), buffer.begin());
}

void MultiConditionProblem::fillInitialParameters(gsl::span<double> buffer) const
{
    if(!startingPoint.empty()) {
        RELEASE_ASSERT(buffer.size() == startingPoint.size(), "");
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
     * TODO: need to have current parameters; need to be supplied to intermediate function;
     * need to from optimizer if in line search or actual step
     */

    // validationProblem->evaluateObjectiveFunction()

    return static_cast<int>(stop);
}


void MultiConditionProblem::setInitialParameters(std::vector<double> const& startingPoint)
{
    this->startingPoint = startingPoint;
}

void MultiConditionProblem::setParametersMin(std::vector<double> const& lowerBounds)
{
    parametersMin = lowerBounds;
}

void MultiConditionProblem::setParametersMax(std::vector<double> const& upperBounds)
{
    parametersMax = upperBounds;
}

std::unique_ptr<OptimizationReporter> MultiConditionProblem::getReporter() const
{

    return std::make_unique<OptimizationReporter>(
                costFun.get(),
                std::make_unique<OptimizationResultWriter>(*resultWriter),
                std::make_unique<Logger>(*logger));
}

std::vector<int> MultiConditionProblem::getTrainingData() const
{
    std::vector<int> dataIndices(dataProvider->getNumberOfConditions());
    std::iota(dataIndices.begin(), dataIndices.end(), 0);
    return dataIndices;
}

MultiConditionDataProvider *MultiConditionProblem::getDataProvider() {
    return dataProvider;
}

OptimizationResultWriter *MultiConditionProblem::getResultWriter() { return resultWriter.get(); }

MultiConditionProblemMultiStartOptimizationProblem::MultiConditionProblemMultiStartOptimizationProblem(MultiConditionDataProviderHDF5 *dp,
        OptimizationOptions options,
        OptimizationResultWriter *resultWriter,
        LoadBalancerMaster *loadBalancer, std::unique_ptr<Logger> logger)
    : dp(dp), options(std::move(options)),
      resultWriter(resultWriter), loadBalancer(loadBalancer),
      logger(std::move(logger))
{}

int MultiConditionProblemMultiStartOptimizationProblem::getNumberOfStarts() const { return options.numStarts; }

bool MultiConditionProblemMultiStartOptimizationProblem::restartOnFailure() const { return options.retryOptimization; }

std::unique_ptr<OptimizationProblem> MultiConditionProblemMultiStartOptimizationProblem::getLocalProblem(
        int multiStartIndex) const {
    // generate new OptimizationProblem with data from dp

    RELEASE_ASSERT(dp != nullptr, "");

    std::unique_ptr<MultiConditionProblem> problem;

    if (resultWriter) {
        problem = std::make_unique<MultiConditionProblem>(
                    dp, loadBalancer,
                    logger->getChild(std::string("o") + std::to_string(multiStartIndex)),
                    std::make_unique<OptimizationResultWriter>(*resultWriter));
        problem->getResultWriter()->setRootPath("/multistarts/" + std::to_string(multiStartIndex));
    } else {
        problem = std::make_unique<MultiConditionProblem>(
                    dp, loadBalancer,
                    logger->getChild(std::string("o") + std::to_string(multiStartIndex)), nullptr);
    }
    problem->setOptimizationOptions(options);
    problem->setInitialParameters(parpe::OptimizationOptions::getStartingPoint(dp->getHdf5FileId(), multiStartIndex));

    if(options.hierarchicalOptimization)
        return std::unique_ptr<OptimizationProblem>(
                    new parpe::HierarchicalOptimizationProblemWrapper(std::move(problem), dp));

    return std::move(problem);
}

void printSimulationResult(Logger *logger, int jobId, amici::ReturnData const* rdata, double timeSeconds) {
    bool with_sensi = rdata->sensi >= amici::SensitivityOrder::first;

    logger->logmessage(LOGLVL_DEBUG, "Result for %d: %g (%d) (%d/%d/%.4fs%c)",
                       jobId, rdata->llh, rdata->status,
                       rdata->numsteps[rdata->numsteps.size() - 1],
                       with_sensi?rdata->numstepsB[0]:0,
                       timeSeconds,
                       with_sensi?'+':'-');


    // check for NaNs, only report first
    if (with_sensi) {
        for (int i = 0; i < rdata->np; ++i) {
            if (std::isnan(rdata->sllh[i])) {
                logger->logmessage(LOGLVL_DEBUG, "Gradient contains NaN at %d", i);
                break;
            }

            if (std::isinf(rdata->sllh[i])) {
                logger->logmessage(LOGLVL_DEBUG, "Gradient contains Inf at %d", i);
                break;
            }
        }
    }
}


void saveSimulation(hid_t file_id, std::string const& pathStr, std::vector<double> const& parameters,
                   double llh, gsl::span<double const> gradient, double timeElapsedInSeconds,
                   gsl::span<double const> states, gsl::span<double const> stateSensi,
                   gsl::span<double const> outputs, int jobId,
                   int status, std::string const& label)
{
    // TODO replace by SimulationResultWriter
    const char *fullGroupPath = pathStr.c_str();

    auto lock = hdf5MutexGetLock();

    hdf5CreateOrExtendAndWriteToDouble2DArray(
        file_id, fullGroupPath, "simulationLogLikelihood", &llh, 1);

    hdf5CreateOrExtendAndWriteToInt2DArray(file_id, fullGroupPath, "jobId",
                                           &jobId, 1);

    if (!gradient.empty()) {
        hdf5CreateOrExtendAndWriteToDouble2DArray(
            file_id, fullGroupPath, "simulationLogLikelihoodGradient", gradient.data(),
            parameters.size());
    } else if(!parameters.empty()) {
        double dummyGradient[parameters.size()];
        std::fill_n(dummyGradient, parameters.size(), NAN);
        hdf5CreateOrExtendAndWriteToDouble2DArray(
            file_id, fullGroupPath, "simulationLogLikelihoodGradient",
            dummyGradient, parameters.size());
    }

    if (!parameters.empty())
        hdf5CreateOrExtendAndWriteToDouble2DArray(
            file_id, fullGroupPath, "simulationParameters", parameters.data(), parameters.size());

    hdf5CreateOrExtendAndWriteToDouble2DArray(file_id, fullGroupPath,
                                              "simulationWallTimeInSec",
                                              &timeElapsedInSeconds, 1);

    if (!states.empty())
        hdf5CreateOrExtendAndWriteToDouble2DArray(
            file_id, fullGroupPath, "simulationStates", states.data(), states.size());

    if (!outputs.empty())
        hdf5CreateOrExtendAndWriteToDouble2DArray(
            file_id, fullGroupPath, "simulationObservables", outputs.data(), outputs.size());

    if (!stateSensi.empty())
        hdf5CreateOrExtendAndWriteToDouble3DArray(
            file_id, fullGroupPath, "simulationStateSensitivities", stateSensi.data(),
            stateSensi.size() / parameters.size(), parameters.size());

    hdf5CreateOrExtendAndWriteToInt2DArray(file_id, fullGroupPath,
                                           "simulationStatus", &status, 1);

    hdf5CreateOrExtendAndWriteToString1DArray(file_id, fullGroupPath,
                                           "simulationLabel", label);

    H5Fflush(file_id, H5F_SCOPE_LOCAL);

}

AmiciSimulationRunner::AmiciResultPackageSimple runAndLogSimulation(
        amici::Solver const &solverTemplate,
        amici::Model &model,
        int conditionIdx,
        int jobId,
        MultiConditionDataProvider const *dataProvider,
        OptimizationResultWriter *resultWriter,
        bool logLineSearch,
        Logger* logger)
{
    // wall time  on worker for current simulation
    WallTimer simulationTimer;

    /* Get ExpData with measurement and fixed parameters. Other model parameters
     * and sensitivity options have been set already */
    auto edata = dataProvider->getExperimentalDataForCondition(conditionIdx);

    /* In case of simulation failure, try rerunning with higher error tolerance
     * for a total of maxNumTrials times */
    constexpr int maxNumTrials = 6; // on failure, rerun simulation
    // Error tolerance relaxation factor upon failure
    constexpr double errorRelaxation = 1e2;
    std::unique_ptr<amici::ReturnData> rdata;

    for(int trial = 1; trial <= maxNumTrials; ++trial) {
        /* It is currently not safe to reuse solver if an exception has
         * occurred,so clone every time */
        auto solver = std::unique_ptr<amici::Solver>(solverTemplate.clone());

        if(trial - 1 == maxNumTrials) {
            logger->logmessage(LOGLVL_ERROR,
                               "Simulation trial %d/%d failed. Giving up.",
                               trial, maxNumTrials);
            break;
        }

        if(rdata) {
            /* something went wrong in the previous simulation. until we have
             * better exception handling, we check those fields to deduce where
             * the error occurred
             */
            bool forwardFailed = std::isnan(rdata->x[rdata->x.size() - 1]);
            bool backwardFailed = std::isnan(rdata->llh);

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
                        LOGLVL_WARNING,
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
            rdata = amici::runAmiciSimulation(*solver, edata.get(), model);
        } catch (std::exception const& e) {
            logger->logmessage(
                        LOGLVL_WARNING, "Error during simulation: %s (%d)",
                        e.what(), rdata->status);
            if(rdata->status == AMICI_SUCCESS)
                // shouldn't happen, but just to be safe
                rdata->status = AMICI_ERROR;
            rdata->invalidateLLH();
            rdata->invalidateSLLH();
        }

        if(rdata->status == AMICI_SUCCESS)
            break;
    }
    double timeSeconds = simulationTimer.getTotal();

    printSimulationResult(logger, jobId, rdata.get(), timeSeconds);

    if (resultWriter && (solverTemplate.getSensitivityOrder()
                         > amici::SensitivityOrder::none || logLineSearch)) {
        saveSimulation(resultWriter->getFileId(), resultWriter->getRootPath(),
                       model.getParameters(), rdata->llh, rdata->sllh,
                       timeSeconds, rdata->x, rdata->sx, rdata->y,
                       jobId, rdata->status, logger->getPrefix());
    }

    return AmiciSimulationRunner::AmiciResultPackageSimple {
        rdata->llh,
                timeSeconds,
                (solverTemplate.getSensitivityOrder()
                 > amici::SensitivityOrder::none)
                ? rdata->sllh : std::vector<double>(),
                rdata->y,
                rdata->status
    };
}

FunctionEvaluationStatus getModelOutputs(
        MultiConditionDataProvider *dataProvider,
        LoadBalancerMaster *loadBalancer,
        int maxSimulationsPerPackage,
        OptimizationResultWriter *resultWriter,
        bool logLineSearch,
        gsl::span<const double> parameters,
        std::vector<std::vector<double> > &modelOutput,
        Logger *logger, double *cpuTime)
{
    int errors = 0;

    std::vector<int> dataIndices(dataProvider->getNumberOfConditions());
    std::iota(dataIndices.begin(), dataIndices.end(), 0);

    modelOutput.resize(dataIndices.size());
    auto parameterVector = std::vector<double>(parameters.begin(),
                                               parameters.end());
    auto jobFinished = [&](JobData *job, int dataIdx) { // jobFinished
        // deserialize
        auto results =
                amici::deserializeFromChar<AmiciSummedGradientFunction::ResultMap> (
                    job->recvBuffer.data(), job->recvBuffer.size());
        job->recvBuffer = std::vector<char>(); // free buffer

        for (auto const& result : results) {
            errors += result.second.status;
            modelOutput[result.first] = result.second.modelOutput;
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
                    [&](std::vector<char> &buffer, int jobId) {
                messageHandler(dataProvider, resultWriter, logLineSearch, buffer, jobId);
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
                    std::vector<char> &buffer, int jobId) {

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
        dataProvider->updateSimulationParameters(
                    conditionIdx,
                    workPackage.optimizationParameters,
                    *model);
        Logger logger(workPackage.logPrefix
                      + "c" + std::to_string(conditionIdx));
        auto result = runAndLogSimulation(
                    *solver, *model, conditionIdx, jobId, dataProvider,
                    resultWriter, logLineSearch, &logger);
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

    setSensitivityOptions(gradient.size());
    fval = 0.0;
    if (gradient.size())
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

FunctionEvaluationStatus AmiciSummedGradientFunction::getModelOutputs(
        gsl::span<const double> parameters,
        std::vector<std::vector<double> > &modelOutput,
        Logger *logger, double *cpuTime) const
{
    return parpe::getModelOutputs(dataProvider, loadBalancer,
                                  maxSimulationsPerPackage, resultWriter,
                                  logLineSearch, parameters, modelOutput,
                                  logger, cpuTime);
}

std::vector<std::vector<double> > AmiciSummedGradientFunction::getAllSigmas() const {
    // TODO: some could be parameter-dependent
    return dataProvider->getAllSigmas();
}

std::vector<std::vector<double> > AmiciSummedGradientFunction::getAllMeasurements() const {
    return dataProvider->getAllMeasurements();
}

void AmiciSummedGradientFunction::messageHandler(std::vector<char> &buffer, int jobId) const {
    parpe::messageHandler(dataProvider, resultWriter, logLineSearch, buffer,
                          jobId);
}

amici::ParameterScaling AmiciSummedGradientFunction::getParameterScaling(int parameterIndex) const
{
    // parameterIndex is optimization parameter index,
    // not necessarily model parameter index!
    return dataProvider->getParameterScale(parameterIndex);
}

int AmiciSummedGradientFunction::runSimulations(
        gsl::span<const double> optimizationParameters,
        double &nllh, gsl::span<double> objectiveFunctionGradient,
        std::vector<int> dataIndices, Logger *logger, double *cpuTime) const {

    int errors = 0;

    auto parameterVector = std::vector<double>(
                optimizationParameters.begin(),
                optimizationParameters.end());
    double simulationTimeSec = 0.0;

    AmiciSimulationRunner simRunner(
                parameterVector,
                objectiveFunctionGradient.size()
                ? amici::SensitivityOrder::first
                : amici::SensitivityOrder::none,
                dataIndices,
                [&](JobData *job, int /*jobIdx*/) {
        errors += aggregateLikelihood(*job,
                                      nllh,
                                      objectiveFunctionGradient,
                                      simulationTimeSec);
    }, nullptr,  logger?logger->getPrefix():"");

#ifdef PARPE_ENABLE_MPI
    if (loadBalancer && loadBalancer->isRunning()) {
        // When running simulations (without gradient),
        // send more simulations to each worker
        // to reduce communication overhead
        errors += simRunner.runDistributedMemory(
                    loadBalancer,
                    objectiveFunctionGradient.size()
                    ? maxGradientSimulationsPerPackage
                    : maxSimulationsPerPackage);
    } else {
#endif
        // Adjoint sensitivity analysis in Sundials 2.6.2 is not thread-safe, so run sequentially
        bool noMultiThreadingWithAdjoints =
                (solver->getSensitivityMethod() == amici::SensitivityMethod::adjoint)
                && !objectiveFunctionGradient.empty();
        errors += simRunner.runSharedMemory(
                    [&](std::vector<char> &buffer, int jobId) {
                messageHandler(buffer, jobId);
    }, noMultiThreadingWithAdjoints);
#ifdef PARPE_ENABLE_MPI
    }
#endif
    if(cpuTime)
        *cpuTime = simulationTimeSec;

    return errors;
}

int AmiciSummedGradientFunction::aggregateLikelihood(JobData &data, double &negLogLikelihood, gsl::span<double> negLogLikelihoodGradient, double &simulationTimeInS) const
{
    int errors = 0;

    // deserialize
    auto results =
            amici::deserializeFromChar<ResultMap> (
                data.recvBuffer.data(), data.recvBuffer.size());
    data.recvBuffer = std::vector<char>(); // free buffer

    for (auto const& result : results) {
        errors += result.second.status != AMICI_SUCCESS;

        // sum up
        negLogLikelihood -= result.second.llh;
        simulationTimeInS += result.second.simulationTimeSeconds;

        if (negLogLikelihoodGradient.size())
            addSimulationGradientToObjectiveFunctionGradient(
                        result.first, result.second.gradient,
                        negLogLikelihoodGradient);

    }
    return errors;
}

void AmiciSummedGradientFunction::addSimulationGradientToObjectiveFunctionGradient(
        int conditionIdx, gsl::span<const double> simulationGradient,
        gsl::span<double> objectiveFunctionGradient) const {
    dataProvider->mapSimulationToOptimizationVariablesAddMultiply(
                conditionIdx, simulationGradient,
                objectiveFunctionGradient, -1.0);
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
