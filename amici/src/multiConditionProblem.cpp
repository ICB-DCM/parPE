#include "multiConditionProblem.h"

#include "steadystateSimulator.h"
#include <optimizationOptions.h>
#include <logging.h>
#include <misc.h>
#include <hierarchicalOptimization.h>

#include <gsl/gsl-lite.hpp>

#include <cassert>
#include <cstring>
#include <ctime>
#include <numeric>


namespace parpe {

// For debugging:
// skip objective function evaluation completely
//#define NO_OBJ_FUN_EVAL

MultiConditionProblem::MultiConditionProblem(MultiConditionDataProvider *dataProvider)
    : MultiConditionProblem(dataProvider, nullptr) {}

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

void MultiConditionProblem::fillParametersMin(double *buffer) const
{
    std::copy(parametersMin.begin(), parametersMin.end(), buffer);
}

void MultiConditionProblem::fillParametersMax(double *buffer) const
{
    std::copy(parametersMax.begin(), parametersMax.end(), buffer);
}

void MultiConditionProblem::fillInitialParameters(double *buffer) const
{
    if(startingPoint.size())
        std::copy(startingPoint.begin(), startingPoint.end(), buffer);
    else
        OptimizationProblem::fillInitialParameters(buffer);
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

    return stop;
}

//template <typename T>
//JobResultAmiciSimulation AmiciSummedGradientFunction<T>::runAndLogSimulation(amici::Solver &solver,
//                                                                             amici::Model &model,
//                                                                             JobIdentifier path,
//                                                                             int jobId) const {
//    double startTime = MPI_Wtime();

//    // run simulation

//    // update UserData::k for condition-specific variables (no parameter mapping
//    // necessary here, this has been done by master)
//    dataProvider->updateFixedSimulationParameters(path.idxConditions, model);

//    auto edata = dataProvider->getExperimentalDataForCondition(path.idxConditions);

//    auto rdata =  std::make_unique<amici::ReturnData>(solver, &model);

//    try {
//        amici::runAmiciSimulation(solver, edata.get(), rdata.get(), model);
//        *rdata->status = AMICI_SUCCESS;
//        rdata->applyChainRuleFactorToSimulationResults(&model);
//    } catch (std::exception &e) {
//        logmessage(LOGLVL_DEBUG, "Error during simulation: %s (%d)", e.what(), (int)*rdata->status);
//    }

//    double endTime = MPI_Wtime();
//    double timeSeconds = (endTime - startTime);

//    printSimulationResult(path, jobId, rdata.get(), timeSeconds);

//    if (resultWriter && (solver.getSensitivityOrder() > amici::AMICI_SENSI_ORDER_NONE || logLineSearch)) {
//        int iterationsUntilSteadystate = -1;
//        logSimulation(resultWriter->getFileId(), resultWriter->getOptimizationPath(),
//                      model.getParameters(), rdata->llh[0], rdata->sllh,
//                timeSeconds, model.np(), model.nx, rdata->x,
//                rdata->sx, model.ny, rdata->y, jobId,
//                iterationsUntilSteadystate, *rdata->status);
//    }
//    int status = (int)*rdata->status;
//    return JobResultAmiciSimulation(status, std::move(rdata), timeSeconds);
//}


//template <typename T>
//void AmiciSummedGradientFunction<T>::messageHandler(std::vector<char> &buffer,
//                                                    int jobId) const {
//    // unpack
//    JobIdentifier path;
//    auto solver = dataProvider->getSolver();
//    auto model = dataProvider->getModel();
//    JobAmiciSimulation<JobIdentifier> sim(solver.get(), model.get(), &path);
//    sim.deserialize(buffer.data(), buffer.size());

//#if QUEUE_WORKER_H_VERBOSE >= 2
//    int mpiRank;
//    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
//    printf("[%d] Received work. ", mpiRank);
//    fflush(stdout);
//#endif

//    // do work
//    JobResultAmiciSimulation result = runAndLogSimulation(*solver, *model, path, jobId);

//#if QUEUE_WORKER_H_VERBOSE >= 2
//    printf("[%d] Work done. ", mpiRank);
//    fflush(stdout);
//#endif

//    // pack & cleanup
//    delete[] result.rdata->J;
//    result.rdata->J = nullptr;

//    delete[] result.rdata->sigmay;
//    result.rdata->sigmay = nullptr;

//    delete[] result.rdata->ssigmay;
//    result.rdata->ssigmay = nullptr;

//    delete[] result.rdata->sx0;
//    result.rdata->sx0 = nullptr;

//    delete[] result.rdata->x;
//    result.rdata->x = nullptr;

//    delete[] result.rdata->x0;
//    result.rdata->x0 = nullptr;

//    delete[] result.rdata->xdot;
//    result.rdata->xdot = nullptr;

//    // will be required for hierarchical optimization
////    delete[] result.rdata->y;
////    result.rdata->y = nullptr;

//    buffer = amici::serializeToStdVec<JobResultAmiciSimulation>(result);
//}




void MultiConditionProblem::setInitialParameters(std::vector<double> startingPoint)
{
    this->startingPoint = startingPoint;
}

void MultiConditionProblem::setParametersMin(std::vector<double> lowerBounds)
{
    parametersMin = lowerBounds;
}

void MultiConditionProblem::setParametersMax(std::vector<double> upperBounds)
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

    problem->setInitialParameters(options.getStartingPoint(dp->getHdf5FileId(), multiStartIndex));

    if(options.hierarchicalOptimization)
        return std::unique_ptr<OptimizationProblem>(
                    new parpe::HierachicalOptimizationProblemWrapper(std::move(problem), dp));
    else
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
            } else if (std::isinf(rdata->sllh[i])) {
                logmessage(LOGLVL_DEBUG, "Result for %s: contains Inf at %d",
                           pathStrBuf, i);
                break;
            }
        }
    }
}


void logSimulation(hid_t file_id, std::string pathStr, std::vector<double> const& theta, double llh, const double *gradient, double timeElapsedInSeconds, int nTheta, int numStates, double *states, double *stateSensi, int numY, double *y, int jobId, int iterationsUntilSteadystate, int status)
{
    // TODO replace by SimulationResultWriter
    const char *fullGroupPath = pathStr.c_str();

    hdf5CreateOrExtendAndWriteToDouble2DArray(
        file_id, fullGroupPath, "simulationLogLikelihood", &llh, 1);

    hdf5CreateOrExtendAndWriteToInt2DArray(file_id, fullGroupPath, "jobId",
                                           &jobId, 1);

    if (gradient) {
        hdf5CreateOrExtendAndWriteToDouble2DArray(
            file_id, fullGroupPath, "simulationLogLikelihoodGradient", gradient,
            nTheta);
    } else {
        double dummyGradient[nTheta];
        std::fill_n(dummyGradient, nTheta, NAN);
        hdf5CreateOrExtendAndWriteToDouble2DArray(
            file_id, fullGroupPath, "simulationLogLikelihoodGradient",
            dummyGradient, nTheta);
    }

    if (theta.size())
        hdf5CreateOrExtendAndWriteToDouble2DArray(
            file_id, fullGroupPath, "simulationParameters", theta.data(), theta.size());

    hdf5CreateOrExtendAndWriteToDouble2DArray(file_id, fullGroupPath,
                                              "simulationWallTimeInSec",
                                              &timeElapsedInSeconds, 1);

    if (states)
        hdf5CreateOrExtendAndWriteToDouble2DArray(
            file_id, fullGroupPath, "simulationStates", states, numStates);

    if (y)
        hdf5CreateOrExtendAndWriteToDouble2DArray(
            file_id, fullGroupPath, "simulationObservables", y, numY);

    hdf5CreateOrExtendAndWriteToInt2DArray(file_id, fullGroupPath,
                                           "iterationsUntilSteadystate",
                                           &iterationsUntilSteadystate, 1);

    if (stateSensi)
        hdf5CreateOrExtendAndWriteToDouble3DArray(
            file_id, fullGroupPath, "simulationStateSensitivities", stateSensi,
            numStates, nTheta);

    hdf5CreateOrExtendAndWriteToInt2DArray(file_id, fullGroupPath,
                                           "simulationStatus", &status, 1);

    H5Fflush(file_id, H5F_SCOPE_LOCAL);

}

SimulationRunnerSimple::AmiciResultPackageSimple runAndLogSimulation(
        amici::Solver &solver, amici::Model &model, JobIdentifier path, int jobId,
        MultiConditionDataProvider *dataProvider, MultiConditionProblemResultWriter *resultWriter,
        bool logLineSearch)
{
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
                      model.getParameters(), rdata->llh, rdata->sllh.data(),
                      timeSeconds, model.np(), model.nx, rdata->x.data(),
                      rdata->sx.data(), model.ny, rdata->y.data(), jobId,
                      iterationsUntilSteadystate, rdata->status);
    }
    return SimulationRunnerSimple::AmiciResultPackageSimple {
        rdata->llh,
                (solver.getSensitivityOrder() > amici::AMICI_SENSI_ORDER_NONE) ? rdata->sllh : std::vector<double>(),
                rdata->y,
                rdata->status
    };
}

} // namespace parpe
