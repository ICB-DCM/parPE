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

// For debugging:
// skip objective function evaluation completely
//#define NO_OBJ_FUN_EVAL

MultiConditionProblem::MultiConditionProblem(
    MultiConditionDataProvider *dataProvider)
    : MultiConditionProblem(dataProvider, nullptr) {}

MultiConditionProblem::MultiConditionProblem(
    MultiConditionDataProvider *dataProvider, LoadBalancerMaster *loadBalancer)
    : dataProvider(dataProvider), loadBalancer(loadBalancer),
      model(dataProvider->getModel()),
      udata(dataProvider->getUserDataForCondition(0)) {

    if (udata == NULL)
        abort();

    init();
}

void MultiConditionProblem::init() {
    numOptimizationParameters = dataProvider->getNumOptimizationParameters();

    parametersMin = new double[numOptimizationParameters];
    dataProvider->getOptimizationParametersLowerBounds(parametersMin);

    parametersMax = new double[numOptimizationParameters];
    dataProvider->getOptimizationParametersUpperBounds(parametersMax);

    lastOptimizationParameters = new double[numOptimizationParameters];
    lastObjectiveFunctionGradient = new double[numOptimizationParameters];
}

int MultiConditionProblem::evaluateObjectiveFunction(const double *optimiziationVariables, double *objectiveFunctionValue,
    double *objectiveFunctionGradient, double *totalTimeInSec) {
    // run on all data
    int numDataIndices = dataProvider->getNumberOfConditions();
    int dataIndices[numDataIndices];
    for (int i = 0; i < numDataIndices; ++i)
        dataIndices[i] = i;

    return evaluateObjectiveFunction(
        optimiziationVariables, objectiveFunctionValue,
        objectiveFunctionGradient, dataIndices, numDataIndices, totalTimeInSec);
}

int MultiConditionProblem::evaluateObjectiveFunction(
    const double *optimiziationVariables, double *objectiveFunctionValue,
    double *objectiveFunctionGradient, int *dataIndices, int numDataIndices, double *totalTimeInSec) {
#ifdef NO_OBJ_FUN_EVAL
    if (objectiveFunctionGradient)
        for (int i = 0; i < numOptimizationParameters; ++i)
            objectiveFunctionGradient = 0;
    *objectiveFunctionValue = 1;
    return 0;
#endif
    // update parameters that are identical for all simulations
    updateUserDataCommon(optimiziationVariables, objectiveFunctionGradient);

    *objectiveFunctionValue = 0;

    if (objectiveFunctionGradient)
        zeros(objectiveFunctionGradient, numOptimizationParameters);

    int errors =
        runSimulations(optimiziationVariables, objectiveFunctionValue,
                       objectiveFunctionGradient, totalTimeInSec,
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

    static double startTime = 0;
    // Wall time on master. NOTE: This also includes waiting time for the job
    // being sent to workers.
    double duration = startTime ? (getTime() - startTime) : 0;

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
        ((MultiConditionProblemResultWriter *)resultWriter)->setJobId(path);
        resultWriter->logLocalOptimizerIteration(
            iter_count, lastOptimizationParameters, numOptimizationParameters,
            obj_value, lastObjectiveFunctionGradient, duration, alg_mod, inf_pr,
            inf_du, mu, d_norm, regularization_size, alpha_du, alpha_pr,
            ls_trials);
    }

    // save start time of the following iteration
    startTime = getTime();

    stop = stop || earlyStopping();

    return stop;
}

void MultiConditionProblem::logObjectiveFunctionEvaluation(
    const double *parameters, double objectiveFunctionValue,
    const double *objectiveFunctionGradient, int numFunctionCalls,
    double timeElapsed) {
    if (resultWriter)
        resultWriter->logLocalOptimizerObjectiveFunctionEvaluation(
            parameters, numOptimizationParameters, objectiveFunctionValue,
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
                                                numOptimizationParameters,
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
     */

    return stop;
}

JobResultAmiciSimulation MultiConditionProblem::runAndLogSimulation(UserData *udata,
                                                       JobIdentifier path,
                                                       int jobId) {
    double startTime = MPI_Wtime();

    Model *model = dataProvider->getModel();

    // run simulation
    int iterationsUntilSteadystate = -1;

    // update UserData::k for condition-specific variables (no parameter mapping
    // necessary here, this has been done by master)
    dataProvider->updateFixedSimulationParameters(path.idxConditions, udata);

    ExpData *edata = dataProvider->getExperimentalDataForCondition(
        path.idxConditions, udata);

    ReturnData *rdata = getSimulationResults(model, udata, edata);
    delete edata;

    double endTime = MPI_Wtime();
    double timeSeconds = (endTime - startTime);

    char pathStrBuf[100];
    path.sprint(pathStrBuf);
    logmessage(LOGLVL_DEBUG, "Result for %s (%d): %g (%d) (%.4fs)", pathStrBuf,
               jobId, rdata->llh[0], (int)*rdata->status, timeSeconds);

    // check for NaNs
    if (udata->sensi >= AMICI_SENSI_ORDER_FIRST) {
        for (int i = 0; i < model->np; ++i)
            if (std::isnan(rdata->sllh[i]))
                logmessage(LOGLVL_DEBUG, "Result for %s: contains NaN at %d",
                           pathStrBuf, i);
            else if (std::isinf(rdata->sllh[i]))
                logmessage(LOGLVL_DEBUG, "Result for %s: contains Inf at %d",
                           pathStrBuf, i);
    }

    if (resultWriter)
        resultWriter->logSimulation(path, udata->p, rdata->llh[0], rdata->sllh,
                                    timeSeconds, model->np, model->nx, rdata->x,
                                    rdata->sx, model->ny, rdata->y, jobId,
                                    iterationsUntilSteadystate, *rdata->status);

    return JobResultAmiciSimulation((int)*rdata->status, rdata, timeSeconds);
}

MultiConditionProblem::~MultiConditionProblem() {
    delete udata;

    delete[] initialParameters;
    delete[] parametersMax;
    delete[] parametersMin;
    delete[] lastObjectiveFunctionGradient;
    delete[] lastOptimizationParameters;
}

void MultiConditionProblem::messageHandler(char **buffer, int *msgSize,
                                           int jobId) {
    // unpack
    UserData *udata = dataProvider->getUserData();
    JobIdentifier path;
    JobAmiciSimulation::toUserData(*buffer, udata, &path);
    delete[] * buffer;

#if QUEUE_WORKER_H_VERBOSE >= 2
    int mpiRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    printf("[%d] Received work. ", mpiRank);
    fflush(stdout);
#endif

    // do work
    JobResultAmiciSimulation result = runAndLogSimulation(udata, path, jobId);

#if QUEUE_WORKER_H_VERBOSE >= 2
    printf("[%d] Work done. ", mpiRank);
    fflush(stdout);
#endif

    // pack & cleanup
    *buffer = serializeToChar<JobResultAmiciSimulation>(&result, msgSize);
    delete result.rdata;
    delete udata;
}

double MultiConditionProblem::getTime() const {
    std::time_t result = std::time(nullptr);
    return result;

    // return MPI_Wtime();
}

double *MultiConditionProblem::getInitialParameters(int multiStartIndex) const {
    return OptimizationOptions::getStartingPoint(dataProvider->getHdf5FileId(),
                                                 multiStartIndex);
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
                                          double *timeInSec,
                                          int *dataIndices,
                                          int numDataIndices) {

    JobIdentifier path = this->path;

    SimulationRunner simRunner(
        [&](int simulationIdx) {
            // extract parameters for simulation of current condition, instead
            // of
            // sending whole  optimization parameter vector to worker
            dataProvider->updateConditionSpecificSimulationParameters(
                dataIndices[simulationIdx], optimizationVariables, udata);
            return udata;
        },
        [&](int simulationIdx) {
            path.idxConditions = dataIndices[simulationIdx];
            return path;
        },
        [&](JobData *jobs, int numJobsTotal) {
            return aggregateLikelihood(jobs, logLikelihood,
                                       objectiveFunctionGradient, timeInSec,
                                       dataIndices,
                                       numDataIndices);
        });

    int errors;

    if (loadBalancer) {
        errors = simRunner.run(
            numDataIndices,
            JobAmiciSimulation::getLength(model->np, sizeof(JobIdentifier)),
            loadBalancer);
    } else {
        errors = simRunner.runSerial(
            numDataIndices,
            JobAmiciSimulation::getLength(model->np, sizeof(JobIdentifier)),
            [&](char **buffer, int *msgSize, int jobId) {
                messageHandler(buffer, msgSize, jobId);
            });
    }

    return errors;
}

int MultiConditionProblem::aggregateLikelihood(
    JobData *data, double *logLikelihood, double *objectiveFunctionGradient, double *timeInSec,
    int *dataIndices, int numDataIndices) {
    int errors = 0;

    // temporary variables for deserialization of simulation results

    for (int simulationIdx = 0; simulationIdx < numDataIndices;
         ++simulationIdx) {
        // deserialize
        JobResultAmiciSimulation result =
                deserializeFromChar<JobResultAmiciSimulation>(
                    data[simulationIdx].recvBuffer,
                    data[simulationIdx].lenRecvBuffer);
        delete[] data[simulationIdx].recvBuffer;
        errors += result.status;
        // sum up
        *logLikelihood -= *result.rdata->llh;
        if(timeInSec)
            *timeInSec += result.simulationTimeInSec;

        if (objectiveFunctionGradient)
            addSimulationGradientToObjectiveFunctionGradient(
                dataIndices[simulationIdx], result.rdata->sllh, objectiveFunctionGradient,
                dataProvider->getNumCommonParameters());
        delete result.rdata;
    }

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
        udata->sensi = AMICI_SENSI_ORDER_FIRST;
        udata->sensi_meth = AMICI_SENSI_ASA;
    } else {
        udata->sensi = AMICI_SENSI_ORDER_NONE;
        udata->sensi_meth = AMICI_SENSI_NONE;
    }
}

void MultiConditionProblem::storeCurrentFunctionEvaluation(
    const double *optimizationParameters, double objectiveFunctionValue,
    const double *objectiveFunctionGradient) {
    lastObjectiveFunctionValue = objectiveFunctionValue;
    memcpy(lastOptimizationParameters, optimizationParameters,
           sizeof(double) * numOptimizationParameters);
    if (objectiveFunctionGradient)
        memcpy(lastObjectiveFunctionGradient, objectiveFunctionGradient,
               sizeof(double) * numOptimizationParameters);
}

MultiConditionDataProvider *MultiConditionProblem::getDataProvider() {
    return dataProvider;
}

OptimizationProblem *
MultiConditionProblemMultiStartOptimization::getLocalProblemImpl(
    int multiStartIndex) {
    // generate new OptimizationProblem with data from dp

    assert(dp != nullptr);
    assert(dp->getModel() != nullptr);

    // Create parallel or serial problem depending on how many processes we are
    // running
    int mpiCommSize = 1;
    int mpiInitialized = 0;
    MPI_Initialized(&mpiInitialized);
    if (mpiInitialized)
        MPI_Comm_size(MPI_COMM_WORLD, &mpiCommSize);

    MultiConditionProblem *problem =
        new MultiConditionProblem(dp, loadBalancer);

    problem->optimizationOptions = new OptimizationOptions(*options);

    if (resultWriter) {
        JobIdentifier id = resultWriter->getJobId();
        id.idxLocalOptimization = multiStartIndex;

        problem->resultWriter =
            new MultiConditionProblemResultWriter(resultWriter->file_id, id);
        problem->path.idxLocalOptimization = multiStartIndex;
    }

    problem->setInitialParameters(
        problem->getInitialParameters(multiStartIndex));

    return problem;
}
