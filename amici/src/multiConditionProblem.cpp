#include "multiConditionProblem.h"

#include "steadystateSimulator.h"
#include <optimizationOptions.h>
#include <logging.h>
#include <misc.h>
#include "hierarchicalOptimization.h"

#include <gsl/gsl-lite.hpp>

#include <cassert>
#include <cstring>
#include <ctime>
#include <numeric>
#include <utility>


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

    costFun = std::make_unique<parpe::SummedGradientFunctionGradientFunctionAdapter<int>>(
                                                                                             std::make_unique<AmiciSummedGradientFunction<int>>(dataProvider, loadBalancer, this->resultWriter.get()),
                                                                                             dataIndices
                                                                                             );
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

MultiConditionProblemMultiStartOptimizationProblem::MultiConditionProblemMultiStartOptimizationProblem(MultiConditionDataProviderHDF5 *dp,
        OptimizationOptions options,
        OptimizationResultWriter *resultWriter,
        LoadBalancerMaster *loadBalancer, std::unique_ptr<Logger> logger)
    : dp(dp), options(std::move(options)),
      resultWriter(resultWriter), loadBalancer(loadBalancer),
      logger(std::move(logger))
{}

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
        problem = std::make_unique<MultiConditionProblem>(dp, loadBalancer,
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
    logger->logmessage(LOGLVL_DEBUG, "Result for %d: %g (%d) (%.4fs%c)",
               jobId, rdata->llh, rdata->status, timeSeconds, rdata->sensi >= amici::SensitivityOrder::first?'+':'-');


    // check for NaNs, only report first
    if (rdata->sensi >= amici::SensitivityOrder::first) {
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
        amici::Solver &solver, amici::Model &model, int conditionIdx, int jobId,
        MultiConditionDataProvider *dataProvider, OptimizationResultWriter *resultWriter,
        bool logLineSearch, Logger* logger)
{
    /* wall time  on worker for current simulation */
    double startTime = MPI_Wtime();

    // run simulation

    // get ExpData with measurement and fixed parameters
    // (other model parameter have been set already, parameter mapping
    // has been done by master)

    auto edata = dataProvider->getExperimentalDataForCondition(conditionIdx);

    // save tolerances to restore later
    auto absTolOrig = solver.getAbsoluteTolerance();
    auto absTolQuadOrig = solver.getAbsoluteToleranceQuadratures();
    auto relTolOrig = solver.getRelativeTolerance();
    auto relTolQuadOrig = solver.getRelativeToleranceQuadratures();

    /* In case of simulation failure, try rerunning with higher error tolerance */
    constexpr int maxNumTrials = 4; // on failure, rerun simulation
    constexpr double errorRelaxation = 1e2;
    std::unique_ptr<amici::ReturnData> rdata;
    for(int trial = 1; trial <= maxNumTrials; ++trial) {
        try {
            rdata = amici::runAmiciSimulation(solver, edata.get(), model);
        } catch (std::exception const& e) {
            logger->logmessage(LOGLVL_DEBUG, "Error during simulation: %s (%d)", e.what(), rdata->status);
            if(rdata->status == AMICI_SUCCESS)
                rdata->status = AMICI_ERROR;
            rdata->invalidateLLH();
        }

        if(rdata->status == AMICI_SUCCESS)
            break;

        if(trial < maxNumTrials) {
            solver.setAbsoluteTolerance(errorRelaxation * solver.getAbsoluteTolerance());
            solver.setAbsoluteToleranceQuadratures(errorRelaxation * solver.getAbsoluteToleranceQuadratures());
            solver.setRelativeTolerance(errorRelaxation * solver.getRelativeTolerance());
            solver.setRelativeToleranceQuadratures(errorRelaxation * solver.getRelativeToleranceQuadratures());

            logger->logmessage(LOGLVL_WARNING, "Error during simulation (try %d/%d), "
                                               "retrying with relaxed error tolerances (*= %g): "
                                               "abs: %g rel: %g quadAbs: %g quadRel: %g ",
                               trial, maxNumTrials, errorRelaxation,
                               solver.getAbsoluteTolerance(), solver.getRelativeTolerance(),
                               solver.getAbsoluteToleranceQuadratures(),
                               solver.getRelativeToleranceQuadratures());
        } else {
            logger->logmessage(LOGLVL_ERROR, "Simulation trial %d/%d failed. Giving up.", trial, maxNumTrials);
        }
    }

    // reset tolerances (solver might be reused)
    solver.setAbsoluteTolerance(absTolOrig);
    solver.setAbsoluteToleranceQuadratures(absTolQuadOrig);
    solver.setRelativeTolerance(relTolOrig);
    solver.setRelativeToleranceQuadratures(relTolQuadOrig);

    double endTime = MPI_Wtime();
    double timeSeconds = (endTime - startTime);

    printSimulationResult(logger, jobId, rdata.get(), timeSeconds);

    if (resultWriter && (solver.getSensitivityOrder() > amici::SensitivityOrder::none || logLineSearch)) {
        saveSimulation(resultWriter->getFileId(), resultWriter->getRootPath(),
                      model.getParameters(), rdata->llh, rdata->sllh,
                      timeSeconds, rdata->x, rdata->sx, rdata->y,
                      jobId, rdata->status, logger->getPrefix());
    }

    return AmiciSimulationRunner::AmiciResultPackageSimple {
        rdata->llh,
                timeSeconds,
                (solver.getSensitivityOrder() > amici::SensitivityOrder::none) ? rdata->sllh : std::vector<double>(),
                rdata->y,
                rdata->status
    };
}

} // namespace parpe
