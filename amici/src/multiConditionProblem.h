#ifndef PARPE_AMICI_MULTI_CONDITION_PROBLEM_H
#define PARPE_AMICI_MULTI_CONDITION_PROBLEM_H

#include <amici/serialization.h>
#include <boost/serialization/map.hpp>

#include "MultiConditionDataProvider.h"
#include "multiConditionProblemResultWriter.h"
#include <multiStartOptimization.h>
#include <optimizationProblem.h>
#include <LoadBalancerMaster.h>
#include <LoadBalancerWorker.h>
#include "simulationRunner.h"

#include <amici/amici.h>

#include <gsl/gsl-lite.hpp>

#include <memory>
#include <cstdlib>

/** @file Interfaces between AMICI model and parPE optimization problem */

namespace parpe {

class MultiConditionDataProvider;
class MultiConditionProblemResultWriter;

SimulationRunnerSimple::AmiciResultPackageSimple  runAndLogSimulation(
        amici::Solver &solver,
        amici::Model &model,
        JobIdentifier path,
        int jobId,
        MultiConditionDataProvider *dataProvider,
        MultiConditionProblemResultWriter *resultWriter,
        bool logLineSearch);


/**
 * @brief The AmiciSummedGradientFunction class represents a cost function
 * based on simulations of an AMICI model for different datasets
 */

template <typename T>
class AmiciSummedGradientFunction : public SummedGradientFunction<T> {

public:
    /**
     * @brief AmiciSummedGradientFunction
     * @param dataProvider Provides data and settings for AMICI simulations
     * @param loadBalancer LoadBalancerMaster for shared memory parallelism, or nullptr
     * @param resultWriter
     */
    AmiciSummedGradientFunction(MultiConditionDataProvider *dataProvider,
                                LoadBalancerMaster *loadBalancer,
                                MultiConditionProblemResultWriter *resultWriter)
        : dataProvider(dataProvider),
          loadBalancer(loadBalancer),
          model(dataProvider->getModel()),
          solver(dataProvider->getSolver()),
          solverOriginal(solver->clone()),
          resultWriter(resultWriter)
    {
        if(auto env = std::getenv("PARPE_LOG_SIMULATIONS")) {
            logLineSearch = env[0] == '1';
        }

    }

    virtual ~AmiciSummedGradientFunction() = default;

    /**
     * @brief Evaluate cost function on a single condition
     * @param parameters
     * @param dataset
     * @param fval
     * @param gradient
     * @return
     */
    virtual FunctionEvaluationStatus evaluate(
            gsl::span<const double> parameters,
            T dataset,
            double &fval,
            gsl::span<double> gradient) const override
    {
        std::vector<T> datasets(1);
        datasets.at(0) = dataset;
        return evaluate(parameters, datasets, fval, gradient);
    }


    virtual FunctionEvaluationStatus evaluate(
            gsl::span<const double> parameters,
            std::vector<T> datasets,
            double &fval,
            gsl::span<double> gradient) const override
    {
#ifdef NO_OBJ_FUN_EVAL
        if (objectiveFunctionGradient)
            std::fill(objectiveFunctionGradient, objectiveFunctionGradient + numOptimizationParameters_, 0);
        *objectiveFunctionValue = 1;
        return 0;
#endif

        setSensitivityOptions(gradient.size());
        fval = 0.0;
        if (gradient.size())
            std::fill(gradient.begin(), gradient.end(), 0.0);

        int errors = runSimulations(parameters, fval, gradient, datasets);

        if (errors || !std::isfinite(fval)) {
            fval = std::numeric_limits<double>::infinity();
            return functionEvaluationFailure;
        }

        return functionEvaluationSuccess;
    }


    /**
     * @brief Number of optimization parameters
     * @return
     */
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
    virtual FunctionEvaluationStatus getModelOutputs(gsl::span<double const> parameters, std::vector<std::vector<double> > &modelOutput) const {
        int errors = 0;
        //    JobIdentifier path; // TODO = this->path;

        std::vector<int> dataIndices(dataProvider->getNumberOfConditions());
        std::iota(dataIndices.begin(), dataIndices.end(), 0);

        setSensitivityOptions(false);
        modelOutput.resize(dataIndices.size());
        auto parameterVector = std::vector<double>(parameters.begin(), parameters.end());
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


        if (loadBalancer && loadBalancer->isRunning()) {
            errors += simRunner.runDistributedMemory(loadBalancer, maxSimulationsPerPackage);
        } else {
            errors += simRunner.runSharedMemory(
                        [&](std::vector<char> &buffer, int jobId) {
                    messageHandler(buffer, jobId);
        }, true);
        }

        return errors == 0 ? functionEvaluationSuccess : functionEvaluationFailure;
    }

    virtual std::vector<std::vector<double>> getAllSigmas() const {
        return dataProvider->getAllSigmas();
    }

    virtual std::vector<std::vector<double>> getAllMeasurements() const {
        return dataProvider->getAllMeasurements();
    }

    /**
     * @brief Is called by worker processes to run a simulation for the given
     * condition
     * @param model Model for simulation. Sensitivity
     * options and parameters are set.
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
     * @brief Callback function for LoadBalancer
     * @param buffer In/out: message buffer
     * @param jobId: In: Identifier of the job (unique up to INT_MAX)
     */
    virtual void messageHandler(std::vector<char> &buffer, int jobId) const {
        // unpack simulation job data
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
        solver->setSensitivityOrder(sim.sensitivityOrder);

        std::map<int, SimulationRunnerSimple::AmiciResultPackageSimple> results;
        // run simulations for all condition indices
        for(auto conditionIndex: sim.conditionIndices) {
            path.idxConditions = conditionIndex;
            dataProvider->updateSimulationParameters(conditionIndex, sim.optimizationParameters, *model);
            auto result = runAndLogSimulation(*solver, *model, path, jobId);
            results[conditionIndex] = result;
        }

#if QUEUE_WORKER_H_VERBOSE >= 2
        printf("[%d] Work done. ", mpiRank);
        fflush(stdout);
#endif
        // serialize to output buffer
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
    virtual int runSimulations(gsl::span<double const> optimizationParameters,
                               double &nllh,
                               gsl::span<double> objectiveFunctionGradient,
                               std::vector<int> dataIndices) const {

        int errors = 0;

        auto parameterVector = std::vector<double>(optimizationParameters.begin(), optimizationParameters.end());

        SimulationRunnerSimple simRunner(parameterVector,
                                         objectiveFunctionGradient.size()?amici::AMICI_SENSI_ORDER_FIRST:amici::AMICI_SENSI_ORDER_NONE,
                                         dataIndices,
                                         [&](JobData *job, int /*jobIdx*/) {
            double simulationTimeSec = 0.0; // TODO not used
            errors += aggregateLikelihood(*job,
                                          nllh,
                                          objectiveFunctionGradient,
                                          simulationTimeSec);
        }, nullptr);
        if (loadBalancer && loadBalancer->isRunning()) {
            // When running simulations (without gradient), send more simulations to each worker
            // to reduce communication overhead
            errors += simRunner.runDistributedMemory(loadBalancer,
                                                     objectiveFunctionGradient.size() ? maxGradientSimulationsPerPackage : maxSimulationsPerPackage);
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
                            gsl::span<double> negLogLikelihoodGradient, double &simulationTimeInS) const {
        int errors = 0;

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
            simulationTimeInS += result.second.simulationTimeSeconds;

            if (negLogLikelihoodGradient.size())
                addSimulationGradientToObjectiveFunctionGradient(
                            result.first, result.second.gradient,
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

    void addSimulationGradientToObjectiveFunctionGradient(int conditionIdx,
                                                          gsl::span<const double> simulationGradient,
                                                          gsl::span<double> objectiveFunctionGradient) const {
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
    const int maxSimulationsPerPackage = 8;
    const int maxGradientSimulationsPerPackage = 1;
};



/**
 * @brief The MultiConditionProblem class represents an optimization problem based
 * on an MultiConditionGradientFunction (AMICI ODE model)
 */

class MultiConditionProblem : public OptimizationProblem {

  public:
    MultiConditionProblem() = default;

    MultiConditionProblem(MultiConditionDataProvider *dp);

    MultiConditionProblem(MultiConditionDataProvider *dp,
                          LoadBalancerMaster *loadBalancer,
                          std::unique_ptr<MultiConditionProblemResultWriter> resultWriter);

    ~MultiConditionProblem() override = default;

    /**
     * @brief earlyStopping
     * @return stop the optimization run
     */
    virtual int earlyStopping();

    MultiConditionDataProvider *getDataProvider();
    MultiConditionProblemResultWriter *getResultWriter() { return resultWriter.get(); }

    //    virtual std::unique_ptr<double[]> getInitialParameters(int multiStartIndex) const override;

    JobIdentifier path;

    void setInitialParameters(const std::vector<double> &startingPoint);
    void setParametersMin(const std::vector<double> &lowerBounds);
    void setParametersMax(std::vector<double> const& upperBounds);

    void fillParametersMin(gsl::span<double> buffer) const override;
    void fillParametersMax(gsl::span<double> buffer) const override;
    void fillInitialParameters(gsl::span<double> buffer) const override;

    std::unique_ptr<OptimizationReporter> getReporter() const override;

  protected:
    //TODO std::unique_ptr<OptimizationProblem> validationProblem;

    MultiConditionDataProvider *dataProvider = nullptr;

private:
    std::unique_ptr<MultiConditionProblemResultWriter> resultWriter;

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
    MultiConditionProblemMultiStartOptimizationProblem(MultiConditionDataProviderHDF5 *dp,
                                                       OptimizationOptions options,
                                                       MultiConditionProblemResultWriter *resultWriter,
                                                       LoadBalancerMaster *loadBalancer);


    int getNumberOfStarts() const { return options.numStarts; }

    bool restartOnFailure() const { return options.retryOptimization; }

    std::unique_ptr<OptimizationProblem> getLocalProblem(int multiStartIndex) const override;

private:
    MultiConditionDataProviderHDF5 *dp = nullptr;
    OptimizationOptions options;
    MultiConditionProblemResultWriter *resultWriter = nullptr;
    LoadBalancerMaster *loadBalancer = nullptr;
};



void printSimulationResult(JobIdentifier const& path, int jobId, const amici::ReturnData *rdata, double timeSeconds);

void logSimulation(hid_t file_id, const std::string &pathStr, const std::vector<double> &parameters, double llh,
                   gsl::span<const double> gradient, double timeElapsedInSeconds, gsl::span<const double> states,
                   gsl::span<const double> stateSensi, gsl::span<const double> outputs, int jobId,
                   int iterationsUntilSteadystate, int status);


} // namespace parpe

#endif
