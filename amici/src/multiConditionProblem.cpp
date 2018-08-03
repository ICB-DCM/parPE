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
    : MultiConditionProblem(dp, nullptr) {}

MultiConditionProblem::MultiConditionProblem(
        MultiConditionDataProvider *dp, LoadBalancerMaster *loadBalancer)
    :dataProvider(dp)
{
    // run on all data
    std::vector<int> dataIndices(dataProvider->getNumberOfConditions());
    std::iota(dataIndices.begin(), dataIndices.end(), 0);

    costFun = std::make_unique<parpe::SummedGradientFunctionGradientFunctionAdapter<int>>(
                                                                                             std::make_unique<AmiciSummedGradientFunction<int>>(dataProvider, loadBalancer),
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

    return std::make_unique<OptimizationReporter>(costFun.get(),
                std::make_unique<MultiConditionProblemResultWriter>(*resultWriter));
}



MultiConditionDataProvider *MultiConditionProblem::getDataProvider() {
    return dataProvider;
}

MultiConditionProblemMultiStartOptimizationProblem::MultiConditionProblemMultiStartOptimizationProblem(
        MultiConditionDataProviderHDF5 *dp,
        OptimizationOptions options,
        MultiConditionProblemResultWriter *resultWriter,
        LoadBalancerMaster *loadBalancer)
    : dp(dp), options(std::move(options)), resultWriter(resultWriter), loadBalancer(loadBalancer)
{}

std::unique_ptr<OptimizationProblem> MultiConditionProblemMultiStartOptimizationProblem::getLocalProblem(
        int multiStartIndex) const {
    // generate new OptimizationProblem with data from dp

    assert(dp != nullptr);

    std::unique_ptr<MultiConditionProblem> problem = std::make_unique<MultiConditionProblem>(dp, loadBalancer);

    problem->setOptimizationOptions(options);

    if (resultWriter) {
        JobIdentifier id = resultWriter->getJobId();
        id.idxLocalOptimization = multiStartIndex;

        problem->resultWriter = std::make_unique<MultiConditionProblemResultWriter>(*resultWriter);
        problem->resultWriter->setJobId(id);
        problem->path.idxLocalOptimization = multiStartIndex;
    }

    problem->setInitialParameters(parpe::OptimizationOptions::getStartingPoint(dp->getHdf5FileId(), multiStartIndex));

    if(options.hierarchicalOptimization)
        return std::unique_ptr<OptimizationProblem>(
                    new parpe::HierachicalOptimizationProblemWrapper(std::move(problem), dp));

    return std::move(problem);
}

void printSimulationResult(const JobIdentifier &path, int jobId, amici::ReturnData const* rdata, double timeSeconds) {
    char pathStrBuf[100];
    path.sprint(pathStrBuf);
    logmessage(LOGLVL_DEBUG, "Result for %s (%d): %g (%d) (%.4fs)", pathStrBuf,
               jobId, rdata->llh, rdata->status, timeSeconds);


    // check for NaNs
    if (rdata->sensi >= amici::AMICI_SENSI_ORDER_FIRST) {
        for (int i = 0; i < rdata->np; ++i) {
            if (std::isnan(rdata->sllh[i])) {
                logmessage(LOGLVL_DEBUG, "Result for %s: contains NaN at %d",
                           pathStrBuf, i);
                break;
            }

            if (std::isinf(rdata->sllh[i])) {
                logmessage(LOGLVL_DEBUG, "Result for %s: contains Inf at %d",
                           pathStrBuf, i);
                break;
            }
        }
    }
}


void logSimulation(hid_t file_id, std::string const& pathStr, std::vector<double> const& parameters,
                   double llh, gsl::span<double const> gradient, double timeElapsedInSeconds,
                   gsl::span<double const> states, gsl::span<double const> stateSensi,
                   gsl::span<double const> outputs, int jobId, int iterationsUntilSteadystate, int status)
{
    // TODO replace by SimulationResultWriter
    const char *fullGroupPath = pathStr.c_str();

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

    hdf5CreateOrExtendAndWriteToInt2DArray(file_id, fullGroupPath,
                                           "iterationsUntilSteadystate",
                                           &iterationsUntilSteadystate, 1);

    if (!stateSensi.empty())
        hdf5CreateOrExtendAndWriteToDouble3DArray(
            file_id, fullGroupPath, "simulationStateSensitivities", stateSensi.data(),
            stateSensi.size() / parameters.size(), parameters.size());

    hdf5CreateOrExtendAndWriteToInt2DArray(file_id, fullGroupPath,
                                           "simulationStatus", &status, 1);

    H5Fflush(file_id, H5F_SCOPE_LOCAL);

}

SimulationRunnerSimple::AmiciResultPackageSimple runAndLogSimulation(
        amici::Solver &solver, amici::Model &model, JobIdentifier path, int jobId,
        MultiConditionDataProvider *dataProvider, MultiConditionProblemResultWriter *resultWriter,
        bool logLineSearch)
{
    /* wall time  on worker for current simulation */
    double startTime = MPI_Wtime();

    // run simulation

    // get ExpData with measurement and fixed parameters
    // (other model parameter have been set already, parameter mapping
    // has been done by master)

    auto edata = dataProvider->getExperimentalDataForCondition(path.idxConditions);

    std::unique_ptr<amici::ReturnData> rdata;
    try {
        rdata = amici::runAmiciSimulation(solver, edata.get(), model);
    } catch (std::exception &e) {
        logmessage(LOGLVL_DEBUG, "Error during simulation: %s (%d)", e.what(), rdata->status);
        rdata->invalidateLLH();
    }

    double endTime = MPI_Wtime();
    double timeSeconds = (endTime - startTime);

    printSimulationResult(path, jobId, rdata.get(), timeSeconds);

    if (resultWriter && (solver.getSensitivityOrder() > amici::AMICI_SENSI_ORDER_NONE || logLineSearch)) {
        int iterationsUntilSteadystate = -1;
        logSimulation(resultWriter->getFileId(), resultWriter->getOptimizationPath(),
                      model.getParameters(), rdata->llh, rdata->sllh,
                      timeSeconds, rdata->x, rdata->sx, rdata->y,
                      jobId, iterationsUntilSteadystate, rdata->status);
    }

    return SimulationRunnerSimple::AmiciResultPackageSimple {
        rdata->llh,
                timeSeconds,
                (solver.getSensitivityOrder() > amici::AMICI_SENSI_ORDER_NONE) ? rdata->sllh : std::vector<double>(),
                rdata->y,
                rdata->status
    };
}

} // namespace parpe
