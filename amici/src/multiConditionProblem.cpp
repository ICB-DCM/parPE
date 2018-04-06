#include "multiConditionProblem.h"

#include "simulationWorkerAmici.h"
#include "steadystateSimulator.h"

#include <optimizationOptions.h>
#include <logging.h>
#include <misc.h>
#include <hierarchicalOptimization.h>

#include <amici/model.h>
#include <amici/rdata.h>
#include <amici/serialization.h>

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

template <typename T>
SimulationRunnerSimple::AmiciResultPackageSimple AmiciSummedGradientFunction<T>::runAndLogSimulation(amici::Solver &solver,
                                                                             amici::Model &model,
                                                                             JobIdentifier path,
                                                                             int jobId) const {
    double startTime = MPI_Wtime();

    // run simulation

    // update UserData::k for condition-specific variables (no parameter mapping
    // necessary here, this has been done by master)
    dataProvider->updateFixedSimulationParameters(path.idxConditions, model);

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


template <typename T>
void AmiciSummedGradientFunction<T>::messageHandler(std::vector<char> &buffer,
                                                    int jobId) const {
    // unpack
    JobIdentifier path;
    auto solver = dataProvider->getSolver();
    auto model = dataProvider->getModel();
    auto sim = amici::deserializeFromChar<
            SimulationRunnerSimple::AmiciWorkPackageSimple>(buffer.data(), buffer.size());

#if QUEUE_WORKER_H_VERBOSE >= 2
    int mpiRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    printf("[%d] Received work. ", mpiRank);
    fflush(stdout);
#endif

    std::map<int, SimulationRunnerSimple::AmiciResultPackageSimple> results;
    // do work
    for(auto conditionIndex: sim.conditionIndices) {
        solver->setSensitivityOrder(sim.sensitivityOrder);
        path.idxConditions = conditionIndex;
        dataProvider->updateSimulationParameters(conditionIndex, sim.optimizationParameters.data(), *model);
        // TODO do parameter mapping
        SimulationRunnerSimple::AmiciResultPackageSimple result = runAndLogSimulation(*solver, *model, path, jobId);
        results[conditionIndex] = result;
    }
#if QUEUE_WORKER_H_VERBOSE >= 2
    printf("[%d] Work done. ", mpiRank);
    fflush(stdout);
#endif

    buffer = amici::serializeToStdVec(results);
}

template<typename T>
amici::AMICI_parameter_scaling AmiciSummedGradientFunction<T>::getParameterScaling(int parameterIndex) const
{
    // parameterIndex is optimization parameter index, not necessarily model parameter index!
    return dataProvider->getParameterScale(parameterIndex);
}


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

    return std::make_unique<OptimizationReporter>(
                std::make_unique<MultiConditionProblemResultWriter>(*resultWriter));
}


template<typename T>
void AmiciSummedGradientFunction<T>::updateUserDataCommon(
        const double *simulationParameters,
        const double *objectiveFunctionGradient) const {
    setSensitivityOptions(objectiveFunctionGradient);

    // update common parameters in UserData, cell-line specific ones are updated
    // later
    auto p = model->getParameters();
    std::copy(simulationParameters, simulationParameters + model->np(), p.data());
    model->setParameters(p);
}


template <typename T>
int AmiciSummedGradientFunction<T>::runSimulations(const double *optimizationVariables,
                                                   double &logLikelihood,
                                                   double *objectiveFunctionGradient,
                                                   std::vector<int> dataIndices) const {

    int errors = 0;
    JobIdentifier path; // TODO = this->path;

//    SimulationRunner simRunner(
//                dataIndices.size(),
//                [&](int simulationIdx) {
//        // extract parameters for simulation of current condition, instead
//        // of sending whole  optimization parameter vector to worker
//        auto myModel = std::unique_ptr<amici::Model>(model->clone());
//        dataProvider->updateConditionSpecificSimulationParameters(
//                    dataIndices[simulationIdx], optimizationVariables, *myModel);
//        return std::pair<std::unique_ptr<amici::Model>,std::unique_ptr<amici::Solver>>(std::move(myModel), std::unique_ptr<amici::Solver>(solver->clone()));
//    },
//    [&](int simulationIdx) {
//        path.idxConditions = dataIndices[simulationIdx];
//        return path;
//    },
//    [&](JobData *job, int dataIdx) {
//        double simulationTimeSec = 0.0; // TODO not used
//        errors += aggregateLikelihood(*job,
//                                      logLikelihood,
//                                      objectiveFunctionGradient,
//                                      dataIndices[dataIdx], simulationTimeSec);
//    }, nullptr);

    auto parameterVector = std::vector<double>(optimizationVariables, optimizationVariables + numParameters());
    SimulationRunnerSimple simRunner(parameterVector,
                                     objectiveFunctionGradient?amici::AMICI_SENSI_ORDER_FIRST:amici::AMICI_SENSI_ORDER_NONE,
                                     dataIndices,
                                     [&](JobData *job, int jobIdx) {
        double simulationTimeSec = 0.0; // TODO not used
        errors += aggregateLikelihood(*job,
                                      logLikelihood,
                                      objectiveFunctionGradient,
                                      simulationTimeSec);
    }, nullptr);
    if (loadBalancer && loadBalancer->isRunning()) {
        // TODO 8 per package; but check for lower number worker!!
        errors += simRunner.runDistributedMemory(loadBalancer, objectiveFunctionGradient?1:8);
    } else {
        errors += simRunner.runSharedMemory(
                    [&](std::vector<char> &buffer, int jobId) {
                messageHandler(buffer, jobId);
    }, true);
    }

    return errors;
}

template <typename T>
int AmiciSummedGradientFunction<T>::aggregateLikelihood(JobData &data, double &logLikelihood,
                                                        double *objectiveFunctionGradient, double &simulationTimeInS) const {
    int errors = 0;

//    // deserialize
//    JobResultAmiciSimulation result =
//            amici::deserializeFromChar<JobResultAmiciSimulation>(
//                data.recvBuffer.data(), data.recvBuffer.size());
//    data.recvBuffer = std::vector<char>(); // free buffer
//    errors += result.status;

//    // sum up
//    logLikelihood -= *result.rdata->llh;
//    simulationTimeInS += result.simulationTimeInSec;

//    if (objectiveFunctionGradient)
//        addSimulationGradientToObjectiveFunctionGradient(
//                    dataIdx, result.rdata->sllh, objectiveFunctionGradient,
//                    dataProvider->getNumCommonParameters());

    // deserialize
     auto results =
            amici::deserializeFromChar<
             std::map<int, SimulationRunnerSimple::AmiciResultPackageSimple> > (
                data.recvBuffer.data(), data.recvBuffer.size());
    data.recvBuffer = std::vector<char>(); // free buffer

    for (auto const& result : results) {
        errors += result.second.status;

        // sum up
        logLikelihood -= result.second.llh;
//        simulationTimeInS += result.simulationTimeInSec;

        if (objectiveFunctionGradient)
            addSimulationGradientToObjectiveFunctionGradient(
                        result.first, result.second.gradient.data(),
                        objectiveFunctionGradient);

    }
    return errors;
}


template <typename T>
void AmiciSummedGradientFunction<T>::addSimulationGradientToObjectiveFunctionGradient(
        int conditionIdx, const double *simulationGradient,
        double *objectiveFunctionGradient) const {
    dataProvider->mapSimulationToOptimizationVariablesAddMultiply(
                conditionIdx, simulationGradient, objectiveFunctionGradient, -1.0);
}


template <typename T>
void AmiciSummedGradientFunction<T>::setSensitivityOptions(bool sensiRequired) const {
    // sensitivities requested?
    if (sensiRequired) {
        solver->setSensitivityOrder(solverOriginal->getSensitivityOrder());
        solver->setSensitivityMethod(solverOriginal->getSensitivityMethod());
    } else {
        solver->setSensitivityOrder(amici::AMICI_SENSI_ORDER_NONE);
        solver->setSensitivityMethod(amici::AMICI_SENSI_NONE);
    }
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

template <typename T>
AmiciSummedGradientFunction<T>::AmiciSummedGradientFunction(MultiConditionDataProvider *dataProvider, LoadBalancerMaster *loadBalancer, MultiConditionProblemResultWriter *resultWriter)
    : dataProvider(dataProvider),
      loadBalancer(loadBalancer),
      model(dataProvider->getModel()),
      solver(dataProvider->getSolver()),
      solverOriginal(solver->clone()),
      resultWriter(resultWriter)
{
}

template<typename T>
FunctionEvaluationStatus AmiciSummedGradientFunction<T>::evaluate(const double * const parameters, T dataset, double &fval, double *gradient) const
{
    std::vector<T> datasets(1);
    datasets.at(0) = dataset;
    return evaluate(parameters, datasets, fval, gradient);
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

template<typename T>
FunctionEvaluationStatus AmiciSummedGradientFunction<T>::evaluate(const double * const parameters, std::vector<T> datasets, double &fval, double *gradient) const {
#ifdef NO_OBJ_FUN_EVAL
    if (objectiveFunctionGradient)
        std::fill(objectiveFunctionGradient, objectiveFunctionGradient + numOptimizationParameters_, 0);
    *objectiveFunctionValue = 1;
    return 0;
#endif

    // update parameters that are identical for all simulations
    updateUserDataCommon(parameters, gradient);

    fval = 0;

    if (gradient)
        amici::zeros(gradient, numParameters());

    int errors =
            runSimulations(parameters, fval, gradient, datasets);

    if (errors) {
        fval = INFINITY;
    }


    return errors == 0 ? functionEvaluationSuccess : functionEvaluationFailure;


}

template<typename T>
int AmiciSummedGradientFunction<T>::numParameters() const
{
    return dataProvider->getNumOptimizationParameters();
}

template<typename T>
FunctionEvaluationStatus AmiciSummedGradientFunction<T>::getModelOutputs(const double * const parameters, std::vector<std::vector<double> >& modelOutput) const {
    int errors = 0;
//    JobIdentifier path; // TODO = this->path;

    std::vector<int> dataIndices(dataProvider->getNumberOfConditions());
    std::iota(dataIndices.begin(), dataIndices.end(), 0);
    updateUserDataCommon(parameters, nullptr);
    setSensitivityOptions(false);

    modelOutput.resize(dataIndices.size());


    auto parameterVector = std::vector<double>(parameters, parameters + numParameters());
    SimulationRunnerSimple simRunner(parameterVector,
                                     amici::AMICI_SENSI_ORDER_NONE,
                                     dataIndices,
                                     [&](JobData *job, int dataIdx) {
        // deserialize
        auto results =
                amici::deserializeFromChar<
                std::map<int, SimulationRunnerSimple::AmiciResultPackageSimple> > (
                    job->recvBuffer.data(), job->recvBuffer.size());
        job->recvBuffer = std::vector<char>(); // free buffer

        for (auto const& result : results) {
            errors += result.second.status;
            modelOutput[result.first] = result.second.modelOutput;
        }
    }, nullptr);
    /*
    SimulationRunner simRunner(
                dataIndices.size(),
                [&](int simulationIdx) {
        // extract parameters for simulation of current condition, instead
        // of sending whole  optimization parameter vector to worker
        auto myModel = std::unique_ptr<amici::Model>(model->clone());
        dataProvider->updateConditionSpecificSimulationParameters(
                    dataIndices[simulationIdx], parameters, *myModel);
        return std::pair<std::unique_ptr<amici::Model>,std::unique_ptr<amici::Solver>>(std::move(myModel), std::unique_ptr<amici::Solver>(solver->clone()));
    },
    [&](int simulationIdx) {
        path.idxConditions = dataIndices[simulationIdx];
        return path;
    },
    [&](JobData *job, int dataIdx) {
        // deserialize
        JobResultAmiciSimulation result =
                amici::deserializeFromChar<JobResultAmiciSimulation>(
                    job->recvBuffer.data(), job->recvBuffer.size());
        job->recvBuffer = std::vector<char>(); // free buffer
        errors += result.status;

        modelOutput[dataIdx] = std::vector<double>(result.rdata->y, result.rdata->y
                                                   + (result.rdata->nt * result.rdata->nytrue));

    }, nullptr);*/


    if (loadBalancer && loadBalancer->isRunning()) {
        errors += simRunner.runDistributedMemory(loadBalancer, 8);
    } else {
        errors += simRunner.runSharedMemory(
                    [&](std::vector<char> &buffer, int jobId) {
                messageHandler(buffer, jobId);
    }, true);
    }

    return errors == 0 ? functionEvaluationSuccess : functionEvaluationFailure;
}

template<typename T>
std::vector<std::vector<double> > AmiciSummedGradientFunction<T>::getAllMeasurements() const {
    return dataProvider->getAllMeasurements();
}

} // namespace parpe
