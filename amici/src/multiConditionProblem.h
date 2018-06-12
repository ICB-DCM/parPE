#ifndef PROBLEM_H
#define PROBLEM_H

#include "MultiConditionDataProvider.h"
#include <simulationWorkerAmici.h>
#include <multiConditionProblemResultWriter.h>
#include <multiStartOptimization.h>
#include <optimizationProblem.h>
#include <LoadBalancerMaster.h>
#include <LoadBalancerWorker.h>
#include <SimulationRunner.h>

#include <amici/amici.h>
#include <amici/serialization.h>
#include <boost/serialization/map.hpp>

#include <memory>
#include <cmath> //NAN

/** @file Interfaces between AMICI model and parPE optimization problem */

namespace parpe {

class MultiConditionDataProvider;
class MultiConditionProblemResultWriter;

SimulationRunnerSimple::AmiciResultPackageSimple  runAndLogSimulation(
        amici::Solver &solver, amici::Model &model, JobIdentifier path,
        int jobId, MultiConditionDataProvider *dataProvider, MultiConditionProblemResultWriter *resultWriter,
        bool logLineSearch);

/**
 * @brief The AmiciSummedGradientFunction class represents a cost function based on simulations of an AMICI model for different datasets
 */

template <typename T>
class AmiciSummedGradientFunction : public SummedGradientFunction<T> {
public:

    AmiciSummedGradientFunction(MultiConditionDataProvider *dataProvider,
                                LoadBalancerMaster *loadBalancer,
                                MultiConditionProblemResultWriter *resultWriter = nullptr)
        : dataProvider(dataProvider),
          loadBalancer(loadBalancer),
          model(dataProvider->getModel()),
          solver(dataProvider->getSolver()),
          solverOriginal(solver->clone()),
          resultWriter(resultWriter)
    {
    }

    virtual ~AmiciSummedGradientFunction() = default;

    virtual FunctionEvaluationStatus evaluate(
            const double* const parameters,
            T dataset,
            double &fval,
            double* gradient) const override
    {
        std::vector<T> datasets(1);
        datasets.at(0) = dataset;
        return evaluate(parameters, datasets, fval, gradient);
    }


    virtual FunctionEvaluationStatus evaluate(
            const double* const parameters,
            std::vector<T> datasets,
            double &fval,
            double* gradient) const override {
#ifdef NO_OBJ_FUN_EVAL
        if (objectiveFunctionGradient)
            std::fill(objectiveFunctionGradient, objectiveFunctionGradient + numOptimizationParameters_, 0);
        *objectiveFunctionValue = 1;
        return 0;
#endif

        setSensitivityOptions(gradient);
        fval = 0;
        if (gradient)
            std::fill(gradient, gradient + numParameters(), 0.0);

        int errors = runSimulations(parameters, fval, gradient, datasets);

        if (errors) {
            fval = INFINITY;
        }

        return errors == 0 ? functionEvaluationSuccess : functionEvaluationFailure;
    }

    virtual int numParameters() const override
    {
        return dataProvider->getNumOptimizationParameters();
    }

    /**
     * @brief Run simulations (no gradient) with given parameters and collect model outputs
     * @param parameters Model parameters for simulation
     * @param modelOutput in: some vector reference, will be resized.
     * output: Vector of double vectors containing AMICI ReturnData::y (nt x ny, column-major)
     * @return Simulation status
     */


    virtual std::vector<std::vector<double>> getAllSigmas() const {
        return dataProvider->getAllSigmas();
    }

    virtual FunctionEvaluationStatus getModelOutputs(const double * const parameters, std::vector<std::vector<double> > &modelOutput) const {
        int errors = 0;
        //    JobIdentifier path; // TODO = this->path;

        std::vector<int> dataIndices(dataProvider->getNumberOfConditions());
        std::iota(dataIndices.begin(), dataIndices.end(), 0);

        setSensitivityOptions(false);
        modelOutput.resize(dataIndices.size());
        auto parameterVector = std::vector<double>(parameters, parameters + numParameters());
        SimulationRunnerSimple simRunner(parameterVector,
                                         amici::AMICI_SENSI_ORDER_NONE,
                                         dataIndices,
                                         [&](JobData *job, int dataIdx) { // jobFinished
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
        },
        nullptr /* aggregate */);
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

    virtual std::vector<std::vector<double>> getAllMeasurements() const {
        return dataProvider->getAllMeasurements();
    }

    /**
     * @brief Is called by worker processes to run a simulation for the given
     * condition
     * @param udata UserData for simulation. Model dimensions, sensitivity
     * options and UserData::p is set, others are not.
     * E.g. UserData::k has to be update if applicable.
     * @param dataProvider
     * @param path
     * @param jobId
     * @param resultWriter
     * @return
     */
    // TODO does not belong here
    SimulationRunnerSimple::AmiciResultPackageSimple  runAndLogSimulation(
            amici::Solver &solver, amici::Model &model, JobIdentifier path,
            int jobId) const
    {
        return parpe::runAndLogSimulation(solver, model, path, jobId, dataProvider, resultWriter, logLineSearch);

    }


    /**
     * @brief Callback function for loadbalancer
     * @param buffer In/out: message buffer
     * @param msgSize In/out: size (bytes) of bufferobjFunVal
     * @param jobId: In: Identifier of the job (unique up to INT_MAX)
     */
    virtual void messageHandler(std::vector<char> &buffer, int jobId) const {
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

    virtual amici::AMICI_parameter_scaling getParameterScaling(int parameterIndex) const
    {
        // parameterIndex is optimization parameter index, not necessarily model parameter index!
        return dataProvider->getParameterScale(parameterIndex);
    }

protected:// for testing
    AmiciSummedGradientFunction() = default;

    /**
     * @brief Run AMICI simulations for conditions with the given indices
     * @param optimizationVariables
     * @param logLikelihood
     * @param objectiveFunctionGradient
     * @param dataIndices
     * @param numDataIndices
     * @return Simulation status, != 0 indicates failure
     */
    virtual int runSimulations(const double *optimizationVariables,
                               double &nllh,
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
                                          nllh,
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

    /**
     * @brief Aggregates loglikelihood received from workers.
     * @param data Simulation job result
     * @param negLogLikelihood output argument to which *negative* log likelihood is added
     * @param negLogLikelihoodGradient output argument to which *negative* log likelihood gradient is added
     * @param simulationTimeInS unused
     * @return
     */

    int aggregateLikelihood(JobData &data, double &negLogLikelihood,
                            double *negLogLikelihoodGradient, double &simulationTimeInS) const {
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
            negLogLikelihood -= result.second.llh;
            //        simulationTimeInS += result.simulationTimeInSec;

            if (negLogLikelihoodGradient)
                addSimulationGradientToObjectiveFunctionGradient(
                            result.first, result.second.gradient.data(),
                            negLogLikelihoodGradient);

        }
        return errors;
    }


    /**
     * @brief Aggregates loglikelihood gradient received from workers.
     * @param conditionIdx
     * @param simulationGradient log-likelihood gradient from simulation
     * @param objectiveFunctionGradient output to which *negative* log-likelihood gradient from simulation is added
     */

    void addSimulationGradientToObjectiveFunctionGradient(int conditionIdx, const double *simulationGradient,
                                                          double *objectiveFunctionGradient) const {
        dataProvider->mapSimulationToOptimizationVariablesAddMultiply(
                    conditionIdx, simulationGradient, objectiveFunctionGradient, -1.0);
    }

    void queueSimulation(JobIdentifier path, JobData *d, int *jobDone,
                         pthread_cond_t *jobDoneChangedCondition,
                         pthread_mutex_t *jobDoneChangedMutex,
                         int lenSendBuffer);

    void setSensitivityOptions(bool sensiRequired) const {
        // sensitivities requested?
        if (sensiRequired) {
            solver->setSensitivityOrder(solverOriginal->getSensitivityOrder());
            solver->setSensitivityMethod(solverOriginal->getSensitivityMethod());
        } else {
            solver->setSensitivityOrder(amici::AMICI_SENSI_ORDER_NONE);
            solver->setSensitivityMethod(amici::AMICI_SENSI_NONE);
        }
    }

private:
    // TODO: make owning
    MultiConditionDataProvider *dataProvider = nullptr;
    // Non-owning
    LoadBalancerMaster *loadBalancer = nullptr;
    std::unique_ptr<amici::Model> model;
    std::unique_ptr<amici::Solver> solver;
    std::unique_ptr<amici::Solver> solverOriginal; // for saving sensitivity options which are changed depending on whether gradient is needed
    MultiConditionProblemResultWriter *resultWriter = nullptr; // TODO: owning?
    bool logLineSearch = false;
};



/**
 * @brief The MultiConditionProblem class represents an optimization problem based
 * on an MultiConditionGradientFunction (AMICI ODE model)
 */

class MultiConditionProblem : public OptimizationProblem {

  public:
    MultiConditionProblem() = default;

    MultiConditionProblem(MultiConditionDataProvider *dataProvider);

    MultiConditionProblem(MultiConditionDataProvider *dataProvider,
                          LoadBalancerMaster *loadBalancer);

    ~MultiConditionProblem() = default;


    /**
     * @brief This function is called after each iteration.
     * @return status code, non-zero to abort optimization
     */
//    virtual int intermediateFunction(int alg_mod, int iter_count,
//                                     double obj_value, double inf_pr,
//                                     double inf_du, double mu, double d_norm,
//                                     double regularization_size,
//                                     double alpha_du, double alpha_pr,
//                                     int ls_trials) override;

    /**
     * @brief Called after each cost function evaluation for logging results.
     * @param parameters
     * @param objectiveFunctionValue
     * @param objectiveFunctionGradient
     * @param numFunctionCalls
     * @param timeElapsed
     */
//    virtual void logObjectiveFunctionEvaluation(
//        const double *parameters, double objectiveFunctionValue,
//        const double *objectiveFunctionGradient, int numFunctionCalls,
//        double timeElapsed) override;

    /**
     * @brief Called at the end of an optimization for logging results
     * @param optimalCost
     * @param optimalParameters
     * @param masterTime
     * @param exitStatus
     */
//    virtual void logOptimizerFinished(double optimalCost,
//                                      const double *optimalParameters,
//                                      double masterTime,
//                                      int exitStatus) override;

    /**
     * @brief earlyStopping
     * @return stop the optimization run
     */
    virtual int earlyStopping();
    MultiConditionDataProvider *getDataProvider();
//    virtual std::unique_ptr<double[]> getInitialParameters(int multiStartIndex) const override;

    JobIdentifier path;

    std::unique_ptr<MultiConditionProblemResultWriter> resultWriter;

    void setInitialParameters(std::vector<double> startingPoint);
    void setParametersMin(std::vector<double> lowerBounds);
    void setParametersMax(std::vector<double> upperBounds);

    void fillParametersMin(double *buffer) const override;
    void fillParametersMax(double *buffer) const override;
    void fillInitialParameters(double *buffer) const override;

    std::unique_ptr<OptimizationReporter> getReporter() const;

  protected:


    /**
     * @brief Keep information from last evaluation to avoid recomputation for
     * same parameters
     * @param optimizationParameters
     * @param objectiveFunctionValue
     * @param objectiveFunctionGradient
     */
    void
    storeCurrentFunctionEvaluation(const double *optimizationParameters,
                                   double objectiveFunctionValue,
                                   const double *objectiveFunctionGradient);
    //TODO
    std::unique_ptr<OptimizationProblem> validationProblem;

    MultiConditionDataProvider *dataProvider = nullptr;

private:
    std::vector<double> startingPoint;
    std::vector<double> parametersMin;
    std::vector<double> parametersMax;

};




/**
 * @brief The MultiConditionProblemGeneratorForMultiStart class generates new
 * MultiConditionProblem instances with proper DataProviders for multi-start
 * optimization
 */

class MultiConditionProblemMultiStartOptimizationProblem
    : public MultiStartOptimizationProblem {
  public:
    int getNumberOfStarts() const { return options.numStarts; }

    bool restartOnFailure() const { return options.retryOptimization; }

    std::unique_ptr<OptimizationProblem> getLocalProblem(int multiStartIndex) const override;

    MultiConditionDataProviderHDF5 *dp = nullptr;
    OptimizationOptions options;
    MultiConditionProblemResultWriter *resultWriter = nullptr;
    LoadBalancerMaster *loadBalancer = nullptr;
};



void printSimulationResult(JobIdentifier const& path, int jobId, const amici::ReturnData *rdata, double timeSeconds);

void logSimulation(hid_t file_id, std::string path, const std::vector<double> &theta, double llh,
                   const double *gradient, double timeElapsedInSeconds,
                   int nTheta, int numStates, double *states,
                   double *stateSensi, int numY, double *y, int jobId,
                   int iterationsUntilSteadystate, int status);


} // namespace parpe

#endif
