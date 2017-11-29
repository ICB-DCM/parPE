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
    MultiConditionDataProvider *dataProvider, LoadBalancerMaster *loadBalancer)
    : OptimizationProblem(dataProvider?dataProvider->getNumOptimizationParameters():0),
      dataProvider(dataProvider), loadBalancer(loadBalancer),
      model(dataProvider->getModel()),
      udata(dataProvider->getUserDataForCondition(0)),
      udataOriginal (*udata) {

    if (udata == NULL)
        abort();

    dataProvider->getOptimizationParametersLowerBounds(parametersMin_.data());
    dataProvider->getOptimizationParametersUpperBounds(parametersMax_.data());

    lastOptimizationParameters.resize(numOptimizationParameters_);
    lastObjectiveFunctionGradient.resize(numOptimizationParameters_);

}


int MultiConditionProblem::evaluateObjectiveFunction(const double *optimiziationVariables, double *objectiveFunctionValue,
    double *objectiveFunctionGradient) {
    // run on all data
    int numDataIndices = dataProvider->getNumberOfConditions();
    int dataIndices[numDataIndices];
    std::iota(dataIndices, dataIndices + numDataIndices, 0);

    return evaluateObjectiveFunction(
        optimiziationVariables, objectiveFunctionValue,
        objectiveFunctionGradient, dataIndices, numDataIndices);
}

int MultiConditionProblem::evaluateObjectiveFunction(
    const double *optimiziationVariables, double *objectiveFunctionValue,
    double *objectiveFunctionGradient, int *dataIndices, int numDataIndices) {
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
        amici::zeros(objectiveFunctionGradient, numOptimizationParameters_);

    int errors =
        runSimulations(optimiziationVariables, objectiveFunctionValue,
                       objectiveFunctionGradient,
                       dataIndices, numDataIndices);

    if (errors) {
        printObjectiveFunctionFailureMessage();
        *objectiveFunctionValue = INFINITY;
    }

    storeCurrentFunctionEvaluation(optimiziationVariables,
                                   *objectiveFunctionValue,
                                   objectiveFunctionGradient);

    return errors;
}

int MultiConditionProblem::intermediateFunction(
    int alg_mod, int iter_count, double obj_value, double inf_pr, double inf_du,
    double mu, double d_norm, double regularization_size, double alpha_du,
    double alpha_pr, int ls_trials) {

//    static double startTime = 0;
    // Wall time on master. NOTE: This also includes waiting time for the job
    // being sent to workers.
//    double duration = startTime ? (getTime() - startTime) : 0;

    bool stop = false;

    // update iteration counter for logging
    path.idxLocalOptimizationIteration = iter_count;

    char strBuf[50];
    path.sprint(strBuf);
    //    logmessage(LOGLVL_INFO, "%s: %d %d %e %e %e %e %e %e %e %e %d",
    //    strBuf,
    //               alg_mod, iter_count, obj_value, inf_pr, inf_du,
    //               mu, d_norm, regularization_size, alpha_du, alpha_pr,
    //               ls_trials);

    if (resultWriter) {
        resultWriter->setJobId(path);
        resultWriter->logLocalOptimizerIteration(
            iter_count, lastOptimizationParameters.data(), numOptimizationParameters_,
            obj_value, lastObjectiveFunctionGradient.data(), simulationTimeInS, alg_mod, inf_pr,
            inf_du, mu, d_norm, regularization_size, alpha_du, alpha_pr,
            ls_trials);
    }

    // save start time of the following iteration
//    startTime = getTime();

    stop = stop || earlyStopping();

    return stop;
}

void MultiConditionProblem::logObjectiveFunctionEvaluation(
    const double *parameters, double objectiveFunctionValue,
    const double *objectiveFunctionGradient, int numFunctionCalls,
    double timeElapsed) {
    if (resultWriter)
        resultWriter->logLocalOptimizerObjectiveFunctionEvaluation(
            parameters, numOptimizationParameters_, objectiveFunctionValue,
            objectiveFunctionGradient, numFunctionCalls, timeElapsed);
}

void MultiConditionProblem::logOptimizerFinished(
    double optimalCost, const double *optimalParameters, double masterTime,
    int exitStatus) {

    char strBuf[100];
    path.sprint(strBuf);
    logmessage(LOGLVL_INFO, "%s: Optimizer status %d, final llh: %e, time: %f.",
               strBuf, exitStatus, optimalCost, masterTime);

    if (resultWriter)
        resultWriter->saveLocalOptimizerResults(optimalCost, optimalParameters,
                                                numOptimizationParameters_,
                                                masterTime, exitStatus);
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

JobResultAmiciSimulation MultiConditionProblem::runAndLogSimulation(amici::UserData *udata,
                                                       JobIdentifier path,
                                                       int jobId) {
    double startTime = MPI_Wtime();

    auto model = dataProvider->getModel();

    // run simulation
    int iterationsUntilSteadystate = -1;

    // update UserData::k for condition-specific variables (no parameter mapping
    // necessary here, this has been done by master)
    dataProvider->updateFixedSimulationParameters(path.idxConditions, *udata);

    auto edata = dataProvider->getExperimentalDataForCondition(path.idxConditions, udata);

    auto rdata = std::unique_ptr<amici::ReturnData>(
                amici::getSimulationResults(model, udata, edata.get()));

    double endTime = MPI_Wtime();
    double timeSeconds = (endTime - startTime);

    printSimulationResult(path, jobId, udata, rdata.get(), timeSeconds);

    if (resultWriter)
        resultWriter->logSimulation(path, udata->p, rdata->llh[0], rdata->sllh,
                                    timeSeconds, model->np, model->nx, rdata->x,
                                    rdata->sx, model->ny, rdata->y, jobId,
                                    iterationsUntilSteadystate, *rdata->status);

    int status = (int)*rdata->status;
    return JobResultAmiciSimulation(status, std::move(rdata), timeSeconds);
}


void MultiConditionProblem::messageHandler(std::vector<char> &buffer,
                                           int jobId) {
    // unpack
    auto udata = std::unique_ptr<amici::UserData>(dataProvider->getUserData());
    JobIdentifier path;
    JobAmiciSimulation::toUserData(buffer.data(), udata.get(), &path);

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

std::unique_ptr<double[]> MultiConditionProblem::getInitialParameters(int multiStartIndex) const {
    return std::unique_ptr<double[]>(OptimizationOptions::getStartingPoint(dataProvider->getHdf5FileId(),
                                                 multiStartIndex));
}

void MultiConditionProblem::updateUserDataCommon(
    const double *simulationParameters,
    const double *objectiveFunctionGradient) {
    setSensitivityOptions(objectiveFunctionGradient);

    // update common parameters in UserData, cell-line specific ones are updated
    // later
    memcpy(udata->p, simulationParameters, sizeof(double) * model->np);
}

int MultiConditionProblem::runSimulations(const double *optimizationVariables,
                                          double *logLikelihood,
                                          double *objectiveFunctionGradient,
                                          int *dataIndices,
                                          int numDataIndices) {

    int errors = 0;
    JobIdentifier path = this->path;

    SimulationRunner simRunner(
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
            errors += aggregateLikelihood(*job,
                                       logLikelihood,
                                       objectiveFunctionGradient,
                                       dataIndices[dataIdx]);
        }, nullptr);


    if (loadBalancer && loadBalancer->isRunning()) {
        errors += simRunner.runMPI(
            numDataIndices,
            JobAmiciSimulation::getLength(model->np, sizeof(JobIdentifier)),
            loadBalancer);
    } else {
        errors += simRunner.runSerial(
            numDataIndices,
            JobAmiciSimulation::getLength(model->np, sizeof(JobIdentifier)),
            [&](std::vector<char> &buffer, int jobId) {
                messageHandler(buffer, jobId);
            });
    }

    return errors;
}

int MultiConditionProblem::aggregateLikelihood(JobData &data, double *logLikelihood,
                        double *objectiveFunctionGradient, int dataIdx) {
    int errors = 0;

    // deserialize
    JobResultAmiciSimulation result =
            amici::deserializeFromChar<JobResultAmiciSimulation>(
                data.recvBuffer,
                data.lenRecvBuffer);
    delete[] data.recvBuffer;
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


void MultiConditionProblem::printObjectiveFunctionFailureMessage() {
    char strBuf[100];
    path.sprint(strBuf);
    logmessage(LOGLVL_ERROR, "%s: Objective function evaluation failed!",
               strBuf);
}

void MultiConditionProblem::addSimulationGradientToObjectiveFunctionGradient(
    int conditionIdx, const double *simulationGradient,
    double *objectiveFunctionGradient, int numCommon) {
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

void MultiConditionProblem::
    addSimulationGradientToObjectiveFunctionGradientConditionSpecificParameters(
        const double *simulationGradient, double *objectiveFunctionGradient,
        int numCommon, int numConditionSpecificParams,
        int firstIndexOfCurrentConditionsSpecificOptimizationParameters) {
    // condition specific parameters: map simulation to optimization parameters
    for (int paramIdx = 0; paramIdx < numConditionSpecificParams; ++paramIdx) {
        int idxOpt =
            firstIndexOfCurrentConditionsSpecificOptimizationParameters +
            paramIdx;
        int idxSim = numCommon + paramIdx;
        objectiveFunctionGradient[idxOpt] -= simulationGradient[idxSim];
    }
}


void MultiConditionProblem::setSensitivityOptions(bool sensiRequired) {
    // sensitivities requested?
    if (sensiRequired) {
        udata->sensi = udataOriginal.sensi;
        udata->sensi_meth = udataOriginal.sensi_meth;
    } else {
        udata->sensi = amici::AMICI_SENSI_ORDER_NONE;
        udata->sensi_meth = amici::AMICI_SENSI_NONE;
    }
}

void MultiConditionProblem::storeCurrentFunctionEvaluation(
    const double *optimizationParameters, double objectiveFunctionValue,
    const double *objectiveFunctionGradient) {
    lastObjectiveFunctionValue = objectiveFunctionValue;
    std::copy(optimizationParameters, optimizationParameters + numOptimizationParameters_, lastOptimizationParameters.begin());
    if (objectiveFunctionGradient)
        std::copy(objectiveFunctionGradient, objectiveFunctionGradient + numOptimizationParameters_, lastObjectiveFunctionGradient.begin());
}

MultiConditionDataProvider *MultiConditionProblem::getDataProvider() {
    return dataProvider;
}

std::unique_ptr<OptimizationProblem> MultiConditionProblemMultiStartOptimization::getLocalProblemImpl(
    int multiStartIndex) {
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

    auto startingPoint = problem->getInitialParameters(multiStartIndex);
    problem->setInitialParameters(startingPoint.get());

    return std::move(problem);
}

void printSimulationResult(const JobIdentifier &path, int jobId, amici::UserData const* udata, amici::ReturnData const* rdata, double timeSeconds) {
    char pathStrBuf[100];
    path.sprint(pathStrBuf);
    logmessage(LOGLVL_DEBUG, "Result for %s (%d): %g (%d) (%.4fs)", pathStrBuf,
               jobId, rdata->llh[0], (int)*rdata->status, timeSeconds);


    // check for NaNs
    if (udata->sensi >= amici::AMICI_SENSI_ORDER_FIRST) {
        for (int i = 0; i < udata->np; ++i)
            if (std::isnan(rdata->sllh[i])) {
                logmessage(LOGLVL_DEBUG, "Result for %s: contains NaN at %d",
                           pathStrBuf, i);
                break;
            }
            else if (std::isinf(rdata->sllh[i])) {
                logmessage(LOGLVL_DEBUG, "Result for %s: contains Inf at %d",
                           pathStrBuf, i);
                break;
            }
    }
}

} // namespace parpe
