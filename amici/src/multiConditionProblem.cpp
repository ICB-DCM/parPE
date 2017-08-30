#include "multiConditionProblem.h"
#include "multiConditionProblemResultWriter.h"
#include "optimizationOptions.h"
#include "simulationWorkerAmici.h"
#include "steadystateSimulator.h"
#include <LoadBalancerMaster.h>
#include <amici_interface_cpp.h>
#include <amici_model.h>
#include <cassert>
#include <cstring>
#include <ctime>
#include <logging.h>
#include <misc.h>
#include <rdata.h>
#include <udata.h>

// For debugging:
// skip objective function evaluation completely
//#define NO_OBJ_FUN_EVAL

MultiConditionProblem::MultiConditionProblem() // for testing only
{}

MultiConditionProblem::MultiConditionProblem(
    MultiConditionDataProvider *dataProvider, LoadBalancerMaster *loadBalancer)
    : MultiConditionProblem() {
    this->dataProvider = dataProvider;
    this->loadBalancer = loadBalancer;
    model = dataProvider->getModel();
    udata = dataProvider->getUserDataForCondition(0);

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

int MultiConditionProblem::evaluateObjectiveFunction(
    const double *optimiziationVariables, double *objectiveFunctionValue,
    double *objectiveFunctionGradient) {
    // run on all data
    int numDataIndices = dataProvider->getNumberOfConditions();
    int dataIndices[numDataIndices];
    for (int i = 0; i < numDataIndices; ++i)
        dataIndices[i] = i;

    return evaluateObjectiveFunction(
        optimiziationVariables, objectiveFunctionValue,
        objectiveFunctionGradient, dataIndices, numDataIndices);
}

int MultiConditionProblem::evaluateObjectiveFunction(
    const double *optimiziationVariables, double *objectiveFunctionValue,
    double *objectiveFunctionGradient, int *dataIndices, int numDataIndices) {
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
                       objectiveFunctionGradient, dataIndices, numDataIndices);

    if (errors) {
        printObjectiveFunctionFailureMessage();
        *objectiveFunctionValue = INFINITY;
    }

    lastObjectiveFunctionValue = *objectiveFunctionValue;
    memcpy(lastOptimizationParameters, optimiziationVariables,
           sizeof(double) * numOptimizationParameters);
    if (objectiveFunctionGradient)
        memcpy(lastObjectiveFunctionGradient, objectiveFunctionGradient,
               sizeof(double) * numOptimizationParameters);

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

    path.idxLocalOptimizationIteration = iter_count;

    char strBuf[50];
    sprintJobIdentifier(strBuf, path);
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

    startTime = getTime();

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
    sprintJobIdentifier(strBuf, path);
    logmessage(LOGLVL_INFO, "%s: Optimizer status %d, final llh: %e, time: %f.",
               strBuf, exitStatus, optimalCost, masterTime);

    if (resultWriter)
        resultWriter->saveLocalOptimizerResults(optimalCost, optimalParameters,
                                                numOptimizationParameters,
                                                masterTime, exitStatus);
}

ReturnData *MultiConditionProblem::runAndLogSimulation(
    UserData *udata, MultiConditionDataProvider *dataProvider,
    JobIdentifier path, int jobId,
    MultiConditionProblemResultWriter *resultWriter, int *status) {

    double startTime = MPI_Wtime();

    Model *model = dataProvider->getModel();

    // run simulation
    int iterationsUntilSteadystate = -1;

    // update UserData::k for condition-specific variables (no parameter mapping
    // necessary here, this has been done by master)
    ExpData *edata =
        dataProvider->getExperimentalDataForExperimentAndUpdateFixedParameters(
            path.idxConditions, udata);

    if (edata == NULL) {
        logmessage(LOGLVL_CRITICAL,
                   "Failed to get experiment data. Check data file. Aborting.");
        abort();
    }

    ReturnData *rdata = getSimulationResults(model, udata, edata);

    delete edata;

    double endTime = MPI_Wtime();
    double timeSeconds = (endTime - startTime);

    *status = (int)*rdata->status;

    char pathStrBuf[100];
    sprintJobIdentifier(pathStrBuf, path);
    logmessage(LOGLVL_DEBUG, "Result for %s (%d): %e  (%d) (%.2fs)", pathStrBuf,
               jobId, rdata->llh[0], *status, timeSeconds);

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
                                    iterationsUntilSteadystate, *status);

    return rdata;
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
    Model *model = dataProvider->getModel();

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
    int status = 0;
    ReturnData *rdata = runAndLogSimulation(udata, dataProvider, path, jobId,
                                            resultWriter, &status);

#if QUEUE_WORKER_H_VERBOSE >= 2
    printf("[%d] Work done. ", mpiRank);
    fflush(stdout);
#endif

    // pack & cleanup
    *msgSize = JobResultAmiciSimulation::getLength(model->np);
    *buffer = new char[*msgSize];
    JobResultAmiciSimulation::serialize(rdata, udata, status, *buffer);

    delete rdata;
    delete udata;
}

double MultiConditionProblem::getTime() const { return MPI_Wtime(); }

double *MultiConditionProblem::getInitialParameters(int multiStartIndex) const {
    return OptimizationOptions::getStartingPoint(dataProvider->fileId,
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
                                          int *dataIndices,
                                          int numDataIndices) {
    JobIdentifier path = this->path;

    int numJobsTotal = numDataIndices;
    int numJobsFinished = 0;

    // TODO: allocate and free piecewise or according to max queue length
    JobData *jobs = new JobData[numJobsTotal];

    // TODO: need to send previous steadystate as initial conditions
    int lenSendBuffer =
        JobAmiciSimulation::getLength(model->np, sizeof(JobIdentifier));

    // mutex to wait for simulations to finish
    pthread_cond_t simulationsCond = PTHREAD_COND_INITIALIZER;
    pthread_mutex_t simulationsMutex = PTHREAD_MUTEX_INITIALIZER;

    for (int simulationIdx = 0; simulationIdx < numJobsTotal; ++simulationIdx) {
        // tell worker with condition to work on, for logging and reading proper
        // UserData::k
        path.idxConditions = dataIndices[simulationIdx];

        // extract parameters for simulation of current condition, instead of
        // sending whole  optimization parameter vector to worker
        dataProvider->updateConditionSpecificSimulationParameters(
            dataIndices[simulationIdx], optimizationVariables, udata);

        queueSimulation(path, &jobs[simulationIdx], &numJobsFinished,
                        &simulationsCond, &simulationsMutex, lenSendBuffer);
        // printf("Queued work: "); printDatapath(path);
    }

    // wait for simulations to finish
    pthread_mutex_lock(&simulationsMutex);
    while (numJobsFinished < numJobsTotal) // TODO handle finished simulations
                                           // here, don't wait for all to
                                           // complete; stop early if errors
                                           // occured
        pthread_cond_wait(&simulationsCond, &simulationsMutex);
    pthread_mutex_unlock(&simulationsMutex);
    pthread_mutex_destroy(&simulationsMutex);
    pthread_cond_destroy(&simulationsCond);

    // unpack
    int errors =
        aggregateLikelihood(jobs, logLikelihood, objectiveFunctionGradient,
                            dataIndices, numDataIndices);
    delete[] jobs;

    return errors;
}

int MultiConditionProblem::aggregateLikelihood(
    JobData *data, double *logLikelihood, double *objectiveFunctionGradient,
    int *dataIndices, int numDataIndices) {
    int errors = 0;

    // temporary variables for deserialization of simulation results
    double llhTmp;
    double sllhTmp[model->np];

    for (int simulationIdx = 0; simulationIdx < numDataIndices;
         ++simulationIdx) {
        // deserialize
        errors +=
            unpackSimulationResult(&data[simulationIdx], sllhTmp, &llhTmp);

        // sum up
        *logLikelihood -= llhTmp;

        if (objectiveFunctionGradient)
            addSimulationGradientToObjectiveFunctionGradient(
                dataIndices[simulationIdx], sllhTmp, objectiveFunctionGradient,
                dataProvider->getNumCommonParameters());
    }

    return errors;
}

MultiConditionProblemSerial::MultiConditionProblemSerial() {}

MultiConditionProblemSerial::MultiConditionProblemSerial(
    MultiConditionDataProvider *dataProvider)
    : MultiConditionProblem(dataProvider, NULL) {}

int MultiConditionProblemSerial::runSimulations(
    const double *optimizationVariables, double *logLikelihood,
    double *objectiveFunctionGradient, int *dataIndices, int numDataIndices) {
    JobIdentifier path = this->path;

    int numJobsTotal = numDataIndices;

    // TODO: allocate and free piecewise or according to max queue length
    JobData *data = new JobData[numJobsTotal];

    // TODO: need to send previous steadystate as initial conditions
    int lenSendBuffer =
        JobAmiciSimulation::getLength(model->np, sizeof(JobIdentifier));

    for (int simulationIdx = 0; simulationIdx < numJobsTotal; ++simulationIdx) {
        path.idxConditions = dataIndices[simulationIdx];

        // update condition specific simulation parameters
        dataProvider->updateConditionSpecificSimulationParameters(
            dataIndices[simulationIdx], optimizationVariables, udata);

        data[simulationIdx].lenSendBuffer = lenSendBuffer;
        data[simulationIdx].sendBuffer = new char[lenSendBuffer];

        JobAmiciSimulation work;
        work.data = &path;
        work.lenData = sizeof(path);
        work.sensitivityMethod = udata->sensi_meth;
        work.numSimulationParameters = model->np;
        work.simulationParameters = udata->p;

        work.serialize(data[simulationIdx].sendBuffer);
        int sizeDummy;
        messageHandler(&data[simulationIdx].sendBuffer, &sizeDummy,
                       simulationIdx);
        data[simulationIdx].recvBuffer = data[simulationIdx].sendBuffer;
    }

    // unpack
    int errors =
        aggregateLikelihood(data, logLikelihood, objectiveFunctionGradient,
                            dataIndices, numDataIndices);

    delete[] data;

    return errors;
}

double MultiConditionProblemSerial::getTime() const {
    std::time_t result = std::time(nullptr);
    return result;
}

void MultiConditionProblem::printObjectiveFunctionFailureMessage() {
    char strBuf[100];
    sprintJobIdentifier(strBuf, path);
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

int MultiConditionProblem::unpackSimulationResult(JobData *d,
                                                  double *sllhBuffer,
                                                  double *llh) {
    JobResultAmiciSimulation result;
    result.sllh = sllhBuffer;

    result.deserialize(d->recvBuffer);
    *llh = result.llh;

    delete[] d->recvBuffer;

    return result.status != 0;
}

void MultiConditionProblem::queueSimulation(
    JobIdentifier path, JobData *d, int *jobDone,
    pthread_cond_t *jobDoneChangedCondition,
    pthread_mutex_t *jobDoneChangedMutex, int lenSendBuffer) {
    *d = JobData(lenSendBuffer, new char[lenSendBuffer], jobDone, jobDoneChangedCondition,
                 jobDoneChangedMutex);

    JobAmiciSimulation work;
    work.data = &path;
    work.lenData = sizeof(path);
    work.sensitivityMethod = udata->sensi_meth;
    work.numSimulationParameters = model->np;
    work.simulationParameters = udata->p;
    work.serialize(d->sendBuffer);

    loadBalancer->queueJob(d);
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

MultiConditionDataProvider *MultiConditionProblem::getDataProvider() {
    return dataProvider;
}

OptimizationProblem *
MultiConditionProblemMultiStartOptimization::getLocalProblemImpl(
    int multiStartIndex) {
    // generate new OptimizationProblem with data from dp

    assert(dp != nullptr);
    assert(dp->model != nullptr);

    // Create parallel or serial problem depending on how many processes we are
    // running
    int mpiCommSize = 1;
    int mpiInitialized = 0;
    MPI_Initialized(&mpiInitialized);
    if (mpiInitialized)
        MPI_Comm_size(MPI_COMM_WORLD, &mpiCommSize);

    MultiConditionProblem *problem;

    if (mpiCommSize == 1)
        problem = new MultiConditionProblemSerial(dp);
    else
        problem = new MultiConditionProblem(dp, loadBalancer);

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
