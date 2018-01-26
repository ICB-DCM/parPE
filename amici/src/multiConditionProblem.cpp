#include "multiConditionProblem.h"
#include "multiConditionProblemResultWriter.h"
#include "optimizationOptions.h"
#include "simulationWorkerAmici.h"
#include "steadystateSimulator.h"
#include <LoadBalancerMaster.h>
#include <SimulationRunner.h>
#include <amici_interface_cpp.h>
#include <amici_model.h>
#include <cassert>
#include <cstring>
#include <ctime>
#include <logging.h>
#include <misc.h>
#include <rdata.h>
#include <udata.h>
#include <amici_serialization.h>
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

JobResultAmiciSimulation MultiConditionGradientFunction::runAndLogSimulation(amici::UserData *udata,
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

    if (resultWriter) {
        int iterationsUntilSteadystate = -1;
        resultWriter->logSimulation(path, udata->p, rdata->llh[0], rdata->sllh,
                timeSeconds, model->np, model->nx, rdata->x,
                rdata->sx, model->ny, rdata->y, jobId,
                iterationsUntilSteadystate, *rdata->status);
    }
    int status = (int)*rdata->status;
    return JobResultAmiciSimulation(status, std::move(rdata), timeSeconds);
}


void MultiConditionGradientFunction::messageHandler(std::vector<char> &buffer,
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


void MultiConditionGradientFunction::updateUserDataCommon(
        const double *simulationParameters,
        const double *objectiveFunctionGradient) const {
    setSensitivityOptions(objectiveFunctionGradient);

    // update common parameters in UserData, cell-line specific ones are updated
    // later
    memcpy(udata->p, simulationParameters, sizeof(double) * model->np);
}

int MultiConditionGradientFunction::runSimulations(const double *optimizationVariables,
                                                   double *logLikelihood,
                                                   double *objectiveFunctionGradient,
                                                   int *dataIndices,
                                                   int numDataIndices) const {

    int errors = 0;
    JobIdentifier path; // TODO = this->path;

    SimulationRunner simRunner(
                numDataIndices,
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

int MultiConditionGradientFunction::aggregateLikelihood(JobData &data, double *logLikelihood,
                                                        double *objectiveFunctionGradient, int dataIdx, double &simulationTimeInS) const {
    int errors = 0;

    // deserialize
    JobResultAmiciSimulation result =
            amici::deserializeFromChar<JobResultAmiciSimulation>(
                data.recvBuffer.data(), data.recvBuffer.size());
    data.recvBuffer = std::vector<char>(); // free buffer
    errors += result.status;

    // sum up
    *logLikelihood -= *result.rdata->llh;
    simulationTimeInS += result.simulationTimeInSec;

    if (objectiveFunctionGradient)
        addSimulationGradientToObjectiveFunctionGradient(
                    dataIdx, result.rdata->sllh, objectiveFunctionGradient,
                    dataProvider->getNumCommonParameters());

    return errors;
}



void MultiConditionGradientFunction::addSimulationGradientToObjectiveFunctionGradient(
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

void MultiConditionGradientFunction::
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


void MultiConditionGradientFunction::setSensitivityOptions(bool sensiRequired) const {
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

MultiConditionGradientFunction::MultiConditionGradientFunction(MultiConditionDataProvider *dataProvider, LoadBalancerMaster *loadBalancer, MultiConditionProblemResultWriter *resultWriter)
    : dataProvider(dataProvider), model(dataProvider->getModel()),
      udata(dataProvider->getUserDataForCondition(0)),
      udataOriginal (*udata), resultWriter(resultWriter)
{
    if (udata == nullptr)
        abort();
}

GradientFunction::FunctionEvaluationStatus MultiConditionGradientFunction::evaluate(const double * const parameters, double &objectiveFunctionValue, double *objectiveFunctionGradient) const
{
    // run on all data
    int numDataIndices = dataProvider->getNumberOfConditions();
    int dataIndices[numDataIndices];
    std::iota(dataIndices, dataIndices + numDataIndices, 0);

    return evaluateObjectiveFunction(
                parameters, &objectiveFunctionValue,
                objectiveFunctionGradient, dataIndices, numDataIndices);
}

GradientFunction::FunctionEvaluationStatus MultiConditionGradientFunction::evaluateObjectiveFunction(
        const double *optimiziationVariables, double *objectiveFunctionValue,
        double *objectiveFunctionGradient, int *dataIndices, int numDataIndices) const {
#ifdef NO_OBJ_FUN_EVAL
    if (objectiveFunctionGradient)
        std::fill(objectiveFunctionGradient, objectiveFunctionGradient + numOptimizationParameters_, 0);
    *objectiveFunctionValue = 1;
    return 0;
#endif
    // update parameters that are identical for all simulations
    updateUserDataCommon(optimiziationVariables, objectiveFunctionGradient);

    *objectiveFunctionValue = 0;

    if (objectiveFunctionGradient)
        amici::zeros(objectiveFunctionGradient, numParameters());

    int errors =
            runSimulations(optimiziationVariables, objectiveFunctionValue,
                           objectiveFunctionGradient,
                           dataIndices, numDataIndices);

    if (errors) {
        *objectiveFunctionValue = INFINITY;
    }


    return errors == 0 ? functionEvaluationSuccess : functionEvaluationFailure;
}


int MultiConditionGradientFunction::numParameters() const
{
    return dataProvider->getNumOptimizationParameters();
}

} // namespace parpe
