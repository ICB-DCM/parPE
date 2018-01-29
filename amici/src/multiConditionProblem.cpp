#include "multiConditionProblem.h"

#include "simulationWorkerAmici.h"
#include "steadystateSimulator.h"

#include <optimizationOptions.h>
#include <SimulationRunner.h>
#include <logging.h>
#include <misc.h>

#include <amici_interface_cpp.h>
#include <amici_model.h>
#include <rdata.h>
#include <udata.h>
#include <amici_serialization.h>

#include <cassert>
#include <cstring>
#include <ctime>
#include <numeric>

namespace parpe {

// For debugging:
// skip objective function evaluation completely
//#define NO_OBJ_FUN_EVAL

MultiConditionProblem::MultiConditionProblem(
        MultiConditionDataProvider *dataProvider)
    : MultiConditionProblem(dataProvider, nullptr) {}

MultiConditionProblem::MultiConditionProblem(
        MultiConditionDataProvider *dp, LoadBalancerMaster *loadBalancer)
    : OptimizationProblem(std::unique_ptr<parpe::MultiConditionGradientFunction>(
                              new parpe::MultiConditionGradientFunction(dp, loadBalancer))),
      dataProvider(dp)
{
}

void MultiConditionProblem::fillParametersMin(double *buffer) const
{
    dataProvider->getOptimizationParametersLowerBounds(buffer);
}

void MultiConditionProblem::fillParametersMax(double *buffer) const
{
    dataProvider->getOptimizationParametersUpperBounds(buffer);
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

template <typename T>
JobResultAmiciSimulation AmiciSummedGradientFunction<T>::runAndLogSimulation(amici::UserData *udata,
                                                                             JobIdentifier path,
                                                                             int jobId) const {
    double startTime = MPI_Wtime();

    auto model = dataProvider->getModel();

    // run simulation

    // update UserData::k for condition-specific variables (no parameter mapping
    // necessary here, this has been done by master)
    dataProvider->updateFixedSimulationParameters(path.idxConditions, *udata);

    auto edata = dataProvider->getExperimentalDataForCondition(path.idxConditions, udata);

    auto rdata = std::unique_ptr<amici::ReturnData>(
                amici::getSimulationResults(model, udata, edata.get()));

    double endTime = MPI_Wtime();
    double timeSeconds = (endTime - startTime);

    printSimulationResult(path, jobId, udata, rdata.get(), timeSeconds);

    if (resultWriter && (udata->sensi > amici::AMICI_SENSI_ORDER_NONE || logLineSearch)) {
        int iterationsUntilSteadystate = -1;
        logSimulation(resultWriter->getFileId(), resultWriter->getOptimizationPath(), udata->p, rdata->llh[0], rdata->sllh,
                timeSeconds, model->np, model->nx, rdata->x,
                rdata->sx, model->ny, rdata->y, jobId,
                iterationsUntilSteadystate, *rdata->status);
    }
    int status = (int)*rdata->status;
    return JobResultAmiciSimulation(status, std::move(rdata), timeSeconds);
}


template <typename T>
void AmiciSummedGradientFunction<T>::messageHandler(std::vector<char> &buffer,
                                                    int jobId) const {
    // unpack
    auto udata = std::unique_ptr<amici::UserData>(dataProvider->getUserData());
    JobIdentifier path;
    JobAmiciSimulation<JobIdentifier> sim(udata.get(), &path);
    sim.deserialize(buffer.data(), buffer.size());

#if QUEUE_WORKER_H_VERBOSE >= 2
    int mpiRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    printf("[%d] Received work. ", mpiRank);
    fflush(stdout);
#endif

    // do work
    JobResultAmiciSimulation result = runAndLogSimulation(udata.get(), path, jobId);

#if QUEUE_WORKER_H_VERBOSE >= 2
    printf("[%d] Work done. ", mpiRank);
    fflush(stdout);
#endif

    // pack & cleanup
    delete[] result.rdata->J;
    result.rdata->J = nullptr;

    delete[] result.rdata->sigmay;
    result.rdata->sigmay = nullptr;

    delete[] result.rdata->ssigmay;
    result.rdata->ssigmay = nullptr;

    delete[] result.rdata->sx0;
    result.rdata->sx0 = nullptr;

    delete[] result.rdata->x;
    result.rdata->x = nullptr;

    delete[] result.rdata->x0;
    result.rdata->x0 = nullptr;

    delete[] result.rdata->xdot;
    result.rdata->xdot = nullptr;

    delete[] result.rdata->y;
    result.rdata->y = nullptr;

    buffer = amici::serializeToStdVec<JobResultAmiciSimulation>(&result);
}

double MultiConditionProblem::getTime() const {
    std::time_t result = std::time(nullptr);
    return result;

    // return MPI_Wtime();
}

void MultiConditionProblem::setInitialParameters(std::vector<double> startingPoint)
{
    this->startingPoint = startingPoint;
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
    memcpy(udata->p, simulationParameters, sizeof(double) * model->np);
}


template <typename T>
int AmiciSummedGradientFunction<T>::runSimulations(const double *optimizationVariables,
                                                   double &logLikelihood,
                                                   double *objectiveFunctionGradient,
                                                   std::vector<int> dataIndices) const {

    int errors = 0;
    JobIdentifier path; // TODO = this->path;

    SimulationRunner simRunner(
                dataIndices.size(),
                [&](int simulationIdx) {
        // extract parameters for simulation of current condition, instead
        // of sending whole  optimization parameter vector to worker
        dataProvider->updateConditionSpecificSimulationParameters(
                    dataIndices[simulationIdx], optimizationVariables, udata.get());
        return *udata;
    },
    [&](int simulationIdx) {
        path.idxConditions = dataIndices[simulationIdx];
        return path;
    },
    [&](JobData *job, int dataIdx) {
        double simulationTimeSec = 0.0; // TODO not used
        errors += aggregateLikelihood(*job,
                                      logLikelihood,
                                      objectiveFunctionGradient,
                                      dataIndices[dataIdx], simulationTimeSec);
    }, nullptr);


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

template <typename T>
int AmiciSummedGradientFunction<T>::aggregateLikelihood(JobData &data, double &logLikelihood,
                                                        double *objectiveFunctionGradient, int dataIdx, double &simulationTimeInS) const {
    int errors = 0;

    // deserialize
    JobResultAmiciSimulation result =
            amici::deserializeFromChar<JobResultAmiciSimulation>(
                data.recvBuffer.data(), data.recvBuffer.size());
    data.recvBuffer = std::vector<char>(); // free buffer
    errors += result.status;

    // sum up
    logLikelihood -= *result.rdata->llh;
    simulationTimeInS += result.simulationTimeInSec;

    if (objectiveFunctionGradient)
        addSimulationGradientToObjectiveFunctionGradient(
                    dataIdx, result.rdata->sllh, objectiveFunctionGradient,
                    dataProvider->getNumCommonParameters());

    return errors;
}


template <typename T>
void AmiciSummedGradientFunction<T>::addSimulationGradientToObjectiveFunctionGradient(
        int conditionIdx, const double *simulationGradient,
        double *objectiveFunctionGradient, int numCommon) const {
    // global parameters: simply add
    for (int paramIdx = 0; paramIdx < numCommon; ++paramIdx)
        objectiveFunctionGradient[paramIdx] -= simulationGradient[paramIdx];

    // map condition-specific parameters
    addSimulationGradientToObjectiveFunctionGradientConditionSpecificParameters(
                simulationGradient, objectiveFunctionGradient, numCommon,
                dataProvider->getNumConditionSpecificParametersPerSimulation(),
                dataProvider->getIndexOfFirstConditionSpecificOptimizationParameter(
                    conditionIdx));
}

template <typename T>
void AmiciSummedGradientFunction<T>::
addSimulationGradientToObjectiveFunctionGradientConditionSpecificParameters(
        const double *simulationGradient, double *objectiveFunctionGradient,
        int numCommon, int numConditionSpecificParams,
        int firstIndexOfCurrentConditionsSpecificOptimizationParameters) const {
    // condition specific parameters: map simulation to optimization parameters
    for (int paramIdx = 0; paramIdx < numConditionSpecificParams; ++paramIdx) {
        int idxOpt =
                firstIndexOfCurrentConditionsSpecificOptimizationParameters +
                paramIdx;
        int idxSim = numCommon + paramIdx;
        objectiveFunctionGradient[idxOpt] -= simulationGradient[idxSim];
    }
}

template <typename T>
void AmiciSummedGradientFunction<T>::setSensitivityOptions(bool sensiRequired) const {
    // sensitivities requested?
    if (sensiRequired) {
        udata->sensi = udataOriginal.sensi;
        udata->sensi_meth = udataOriginal.sensi_meth;
    } else {
        udata->sensi = amici::AMICI_SENSI_ORDER_NONE;
        udata->sensi_meth = amici::AMICI_SENSI_NONE;
    }
}


MultiConditionDataProvider *MultiConditionProblem::getDataProvider() {
    return dataProvider;
}

std::unique_ptr<OptimizationProblem> MultiConditionProblemMultiStartOptimizationProblem::getLocalProblem(
        int multiStartIndex) const {
    // generate new OptimizationProblem with data from dp

    assert(dp != nullptr);
    assert(dp->getModel() != nullptr);

    std::unique_ptr<MultiConditionProblem> problem = std::make_unique<MultiConditionProblem>(dp, loadBalancer);

    problem->setOptimizationOptions(options);

    if (resultWriter) {
        JobIdentifier id = resultWriter->getJobId();
        id.idxLocalOptimization = multiStartIndex;

        problem->resultWriter = std::make_unique<MultiConditionProblemResultWriter>(*resultWriter);
        problem->resultWriter->setJobId(id);
        problem->path.idxLocalOptimization = multiStartIndex;
    }

    problem->setInitialParameters(options.getStartingPoint(dp->fileId, multiStartIndex));

    return std::move(problem);
}

void printSimulationResult(const JobIdentifier &path, int jobId, amici::UserData const* udata, amici::ReturnData const* rdata, double timeSeconds) {
    char pathStrBuf[100];
    path.sprint(pathStrBuf);
    logmessage(LOGLVL_DEBUG, "Result for %s (%d): %g (%d) (%.4fs)", pathStrBuf,
               jobId, rdata->llh[0], (int)*rdata->status, timeSeconds);


    // check for NaNs
    if (udata->sensi >= amici::AMICI_SENSI_ORDER_FIRST) {
        for (int i = 0; i < udata->np; ++i) {
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
    : dataProvider(dataProvider), model(dataProvider->getModel()),
      udata(dataProvider->getUserDataForCondition(0)),
      udataOriginal (*udata), resultWriter(resultWriter)
{
    if (udata == nullptr)
        abort();
}

template<typename T>
FunctionEvaluationStatus AmiciSummedGradientFunction<T>::evaluate(const double * const parameters, T dataset, double &fval, double *gradient) const
{
    std::vector<T> datasets(1);
    datasets.at(0) = dataset;
    return evaluate(parameters, datasets, fval, gradient);
}

MultiConditionGradientFunction::MultiConditionGradientFunction(MultiConditionDataProvider *dataProvider,
                                                               LoadBalancerMaster *loadBalancer, MultiConditionProblemResultWriter *resultWriter)
    : summedGradFun(std::make_unique<AmiciSummedGradientFunction<int>>(dataProvider, loadBalancer, resultWriter)),
      numConditions(dataProvider->getNumberOfConditions())
{

}

FunctionEvaluationStatus MultiConditionGradientFunction::evaluate(const double * const parameters, double &objectiveFunctionValue, double *objectiveFunctionGradient) const
{
    // run on all data
    std::vector<int> dataIndices(numConditions);
    std::iota(dataIndices.begin(), dataIndices.end(), 0);

    return summedGradFun->evaluate(parameters, dataIndices, objectiveFunctionValue,
                objectiveFunctionGradient);
}


int MultiConditionGradientFunction::numParameters() const
{
    return summedGradFun->numParameters();
}

void logSimulation(hid_t file_id, std::string pathStr, const double *theta, double llh, const double *gradient, double timeElapsedInSeconds, int nTheta, int numStates, double *states, double *stateSensi, int numY, double *y, int jobId, int iterationsUntilSteadystate, int status)
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

    if (theta)
        hdf5CreateOrExtendAndWriteToDouble2DArray(
            file_id, fullGroupPath, "simulationParameters", theta, nTheta);

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

} // namespace parpe
