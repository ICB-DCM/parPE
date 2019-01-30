#ifndef PARPE_AMICI_MULTI_CONDITION_PROBLEM_H
#define PARPE_AMICI_MULTI_CONDITION_PROBLEM_H

#include "parpeConfig.h"

#include <amici/serialization.h>
#include <boost/serialization/map.hpp>

#include "multiConditionDataProvider.h"
#include <multiStartOptimization.h>
#include <optimizationProblem.h>
#include "amiciSimulationRunner.h"
#include <minibatchOptimization.h>

#ifdef PARPE_ENABLE_MPI
#include <loadBalancerWorker.h>
#endif
#include <loadBalancerMaster.h>

#include <amici/amici.h>

#include <gsl/gsl-lite.hpp>

#include <memory>
#include <cstdlib>

/** @file Interfaces between AMICI model and parPE optimization problem */

namespace parpe {

class MultiConditionDataProvider;

/**
 * @brief Run AMICI simulation for the given condition, save and return results
 * @param solver
 * @param model
 * @param conditionIdx
 * @param jobId
 * @param dataProvider
 * @param resultWriter
 * @param logLineSearch
 * @param logger
 * @return Simulation results
 */
AmiciSimulationRunner::AmiciResultPackageSimple runAndLogSimulation(
        amici::Solver &solver,
        amici::Model &model,
        int conditionIdx,
        int jobId,
        MultiConditionDataProvider *dataProvider,
        OptimizationResultWriter *resultWriter,
        bool logLineSearch,
        Logger *logger);


/**
 * @brief The AmiciSummedGradientFunction class represents a cost function
 * based on simulations of an AMICI model for different datasets
 */

template <typename T>
class AmiciSummedGradientFunction : public SummedGradientFunction<T> {

public:
    using WorkPackage = AmiciSimulationRunner::AmiciWorkPackageSimple;
    using ResultPackage = AmiciSimulationRunner::AmiciResultPackageSimple;
    using ResultMap = std::map<int, ResultPackage>;

    /**
     * @brief AmiciSummedGradientFunction
     * @param dataProvider Provides data and settings for AMICI simulations
     * @param loadBalancer LoadBalancerMaster for shared memory parallelism, or
     * nullptr
     * @param resultWriter
     */
    AmiciSummedGradientFunction(
            MultiConditionDataProvider *dataProvider,
            LoadBalancerMaster *loadBalancer,
            OptimizationResultWriter *resultWriter)
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

        if(auto env = std::getenv("PARPE_MAX_SIMULATIONS_PER_PACKAGE")) {
            maxSimulationsPerPackage = std::stoi(env);
        }

        if(auto env =
                std::getenv("PARPE_MAX_GRADIENT_SIMULATIONS_PER_PACKAGE")) {
            maxGradientSimulationsPerPackage = std::stoi(env);
        }
    }

    virtual ~AmiciSummedGradientFunction() = default;

    virtual FunctionEvaluationStatus evaluate(
            gsl::span<const double> parameters,
            T dataset,
            double &fval,
            gsl::span<double> gradient,
            Logger *logger,
            double *cpuTime) const override
    {
        std::vector<T> datasets(1);
        datasets.at(0) = dataset;
        return evaluate(parameters, datasets, fval, gradient, logger, cpuTime);
    }

    virtual FunctionEvaluationStatus evaluate(
            gsl::span<const double> parameters,
            std::vector<T> datasets,
            double &fval,
            gsl::span<double> gradient,
            Logger *logger,
            double *cpuTime) const override
    {
#ifdef NO_OBJ_FUN_EVAL
        if (objectiveFunctionGradient)
            std::fill(objectiveFunctionGradient,
                      objectiveFunctionGradient + numOptimizationParameters_,
                      0);
        *objectiveFunctionValue = 1;
        return 0;
#endif

        setSensitivityOptions(gradient.size());
        fval = 0.0;
        if (gradient.size())
            std::fill(gradient.begin(), gradient.end(), 0.0);

        int errors = runSimulations(parameters, fval, gradient, datasets,
                                    logger, cpuTime);

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
     * @brief Run simulations (no gradient) with given parameters and collect
     * model outputs
     * @param parameters Model parameters for simulation
     * @param modelOutput in: some vector reference, will be resized.
     * output: Vector of double vectors containing AMICI ReturnData::y
     * (nt x ny, column-major)
     * @return Simulation status
     */
    virtual FunctionEvaluationStatus getModelOutputs(
            gsl::span<double const> parameters,
            std::vector<std::vector<double> > &modelOutput,
            Logger *logger,
            double *cpuTime) const
    {
        int errors = 0;

        std::vector<int> dataIndices(dataProvider->getNumberOfConditions());
        std::iota(dataIndices.begin(), dataIndices.end(), 0);

        setSensitivityOptions(false);
        modelOutput.resize(dataIndices.size());
        auto parameterVector = std::vector<double>(parameters.begin(),
                                                   parameters.end());
        auto jobFinished = [&](JobData *job, int dataIdx) { // jobFinished
            // deserialize
            auto results =
                    amici::deserializeFromChar<ResultMap> (
                        job->recvBuffer.data(), job->recvBuffer.size());
            job->recvBuffer = std::vector<char>(); // free buffer

            for (auto const& result : results) {
                errors += result.second.status;
                modelOutput[result.first] = result.second.modelOutput;
            }
        };
        AmiciSimulationRunner simRunner(parameterVector,
                                        amici::SensitivityOrder::none,
                                        dataIndices,
                                        jobFinished,
                                        nullptr /* aggregate */,
                                        logger?logger->getPrefix():"");


#ifdef PARPE_ENABLE_MPI
        if (loadBalancer && loadBalancer->isRunning()) {
            errors += simRunner.runDistributedMemory(loadBalancer,
                                                     maxSimulationsPerPackage);
        } else {
#endif
            errors += simRunner.runSharedMemory(
                        [&](std::vector<char> &buffer, int jobId) {
                    messageHandler(buffer, jobId);
        });
#ifdef PARPE_ENABLE_MPI
        }
#endif
        return errors == 0 ? functionEvaluationSuccess
                           : functionEvaluationFailure;
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
    AmiciSimulationRunner::AmiciResultPackageSimple  runAndLogSimulation(
            amici::Solver &solver, amici::Model &model, int conditionIdx,
            int jobId, Logger* logger) const
    {
        return parpe::runAndLogSimulation(solver, model,
                                          conditionIdx, jobId, dataProvider,
                                          resultWriter, logLineSearch, logger);
    }


    /**
     * @brief Callback function for LoadBalancer
     * @param buffer In/out: message buffer
     * @param jobId: In: Identifier of the job (unique up to INT_MAX)
     */
    virtual void messageHandler(std::vector<char> &buffer, int jobId) const {

#if QUEUE_WORKER_H_VERBOSE >= 2
        int mpiRank;
        MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
        printf("[%d] Received work. ", mpiRank);
        fflush(stdout);
#endif

        auto solver = dataProvider->getSolver();
        auto model = dataProvider->getModel();

        // unpack simulation job data
        auto workPackage = amici::deserializeFromChar<WorkPackage>(
                    buffer.data(), buffer.size());

        solver->setSensitivityOrder(workPackage.sensitivityOrder);

        ResultMap results;
        // run simulations for all condition indices
        for(auto conditionIdx: workPackage.conditionIndices) {
            dataProvider->updateSimulationParameters(
                        conditionIdx,
                        workPackage.optimizationParameters,
                        *model);
            Logger logger(workPackage.logPrefix
                          + "c" + std::to_string(conditionIdx));
            auto result = runAndLogSimulation(
                        *solver, *model, conditionIdx, jobId, &logger);
            results[conditionIdx] = result;
        }

#if QUEUE_WORKER_H_VERBOSE >= 2
        printf("[%d] Work done. ", mpiRank);
        fflush(stdout);
#endif
        // serialize to output buffer
        buffer = amici::serializeToStdVec(results);
    }

    virtual amici::ParameterScaling getParameterScaling(int parameterIndex) const
    {
        // parameterIndex is optimization parameter index,
        // not necessarily model parameter index!
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
                               std::vector<int> dataIndices,
                               Logger *logger,
                               double *cpuTime) const {

        int errors = 0;

        auto parameterVector = std::vector<double>(
                    optimizationParameters.begin(),
                    optimizationParameters.end());
        double simulationTimeSec = 0.0;

        AmiciSimulationRunner simRunner(
                    parameterVector,
                    objectiveFunctionGradient.size()
                      ? amici::SensitivityOrder::first
                      : amici::SensitivityOrder::none,
                    dataIndices,
                    [&](JobData *job, int /*jobIdx*/) {
            errors += aggregateLikelihood(*job,
                                          nllh,
                                          objectiveFunctionGradient,
                                          simulationTimeSec);
        }, nullptr,  logger?logger->getPrefix():"");

#ifdef PARPE_ENABLE_MPI
        if (loadBalancer && loadBalancer->isRunning()) {
            // When running simulations (without gradient),
            // send more simulations to each worker
            // to reduce communication overhead
            errors += simRunner.runDistributedMemory(
                        loadBalancer,
                        objectiveFunctionGradient.size()
                        ? maxGradientSimulationsPerPackage
                        : maxSimulationsPerPackage);
        } else {
#endif
            // Adjoint sensitivity analysis in Sundials 2.6.2 is not thread-safe, so run sequentially
            bool noMultiThreadingWithAdjoints =
                    (solver->getSensitivityMethod() == amici::SensitivityMethod::adjoint)
                    && !objectiveFunctionGradient.empty();
            errors += simRunner.runSharedMemory(
                        [&](std::vector<char> &buffer, int jobId) {
                    messageHandler(buffer, jobId);
        }, noMultiThreadingWithAdjoints);
#ifdef PARPE_ENABLE_MPI
        }
#endif
        if(cpuTime)
            *cpuTime = simulationTimeSec;

        return errors;
    }

    /**
     * @brief Aggregates loglikelihood received from workers.
     * @param data Simulation job result
     * @param negLogLikelihood output argument to which *negative*
     * log likelihood is added
     * @param negLogLikelihoodGradient output argument to which *negative*
     * log likelihood gradient is added
     * @param simulationTimeInS unused
     * @return
     */

    int aggregateLikelihood(JobData &data, double &negLogLikelihood,
                            gsl::span<double> negLogLikelihoodGradient,
                            double &simulationTimeInS) const
    {
        int errors = 0;

        // deserialize
        auto results =
                amici::deserializeFromChar<ResultMap> (
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
     * @param objectiveFunctionGradient output to which *negative*
     * log-likelihood gradient from simulation is added
     */

    void addSimulationGradientToObjectiveFunctionGradient(
            int conditionIdx,
            gsl::span<const double> simulationGradient,
            gsl::span<double> objectiveFunctionGradient) const {
        dataProvider->mapSimulationToOptimizationVariablesAddMultiply(
                    conditionIdx, simulationGradient,
                    objectiveFunctionGradient, -1.0);
    }

    void setSensitivityOptions(bool sensiRequired) const {
        // sensitivities requested?
        if (sensiRequired) {
            solver->setSensitivityOrder(solverOriginal->getSensitivityOrder());
            solver->setSensitivityMethod(solverOriginal->getSensitivityMethod());
        } else {
            solver->setSensitivityOrder(amici::SensitivityOrder::none);
            solver->setSensitivityMethod(amici::SensitivityMethod::none);
        }
    }

private:
    // TODO: make owning
    MultiConditionDataProvider *dataProvider = nullptr;
    // Non-owning
    LoadBalancerMaster *loadBalancer = nullptr;
    std::unique_ptr<amici::Model> model;
    std::unique_ptr<amici::Solver> solver;
    /** For saving sensitivity options which are changed depending on whether
     * gradient is needed */
    std::unique_ptr<amici::Solver> solverOriginal;
    OptimizationResultWriter *resultWriter = nullptr; // TODO: owning?
    bool logLineSearch = false;
    int maxSimulationsPerPackage = 8;
    int maxGradientSimulationsPerPackage = 1;
};



/**
 * @brief The MultiConditionProblem class represents an optimization problem
 * based on an MultiConditionGradientFunction (AMICI ODE model) and
 * MultiConditionDataProvider
 */

class MultiConditionProblem
        : public MinibatchOptimizationProblem<int>
{

  public:
    MultiConditionProblem() = default;

    MultiConditionProblem(MultiConditionDataProvider *dp);

    MultiConditionProblem(
            MultiConditionDataProvider *dp,
            LoadBalancerMaster *loadBalancer,
            std::unique_ptr<Logger> logger,
            std::unique_ptr<OptimizationResultWriter> resultWriter);

    ~MultiConditionProblem() override = default;

    /**
     * @brief earlyStopping
     * @return stop the optimization run
     */
    virtual int earlyStopping();

    MultiConditionDataProvider *getDataProvider();
    OptimizationResultWriter *getResultWriter() { return resultWriter.get(); }

    //    virtual std::unique_ptr<double[]> getInitialParameters(int multiStartIndex) const override;

    void setInitialParameters(const std::vector<double> &startingPoint);
    void setParametersMin(const std::vector<double> &lowerBounds);
    void setParametersMax(std::vector<double> const& upperBounds);

    void fillParametersMin(gsl::span<double> buffer) const override;
    void fillParametersMax(gsl::span<double> buffer) const override;
    void fillInitialParameters(gsl::span<double> buffer) const override;

    std::unique_ptr<OptimizationReporter> getReporter() const override;

    std::vector<int> getTrainingData() const override;

protected:
    //TODO std::unique_ptr<OptimizationProblem> validationProblem;

    MultiConditionDataProvider *dataProvider = nullptr;

private:
    std::unique_ptr<OptimizationResultWriter> resultWriter;

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
    MultiConditionProblemMultiStartOptimizationProblem(
            MultiConditionDataProviderHDF5 *dp,
            OptimizationOptions options,
            OptimizationResultWriter *resultWriter,
            LoadBalancerMaster *loadBalancer,
            std::unique_ptr<Logger> logger);


    int getNumberOfStarts() const override { return options.numStarts; }

    bool restartOnFailure() const override { return options.retryOptimization; }

    std::unique_ptr<OptimizationProblem> getLocalProblem(int multiStartIndex) const override;

private:
    MultiConditionDataProviderHDF5 *dp = nullptr;
    OptimizationOptions options;
    OptimizationResultWriter *resultWriter = nullptr;
    LoadBalancerMaster *loadBalancer = nullptr;
    std::unique_ptr<Logger> logger;
};


void saveSimulation(
        hid_t file_id, const std::string &pathStr,
        const std::vector<double> &parameters, double llh,
        gsl::span<const double> gradient, double timeElapsedInSeconds,
        gsl::span<const double> states,
        gsl::span<const double> stateSensi, gsl::span<const double> outputs,
        int jobId, int status, const std::string &label);

} // namespace parpe

#endif
