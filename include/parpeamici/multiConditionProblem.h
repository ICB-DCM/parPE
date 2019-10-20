#ifndef PARPE_AMICI_MULTI_CONDITION_PROBLEM_H
#define PARPE_AMICI_MULTI_CONDITION_PROBLEM_H

#include <parpecommon/parpeConfig.h>
#include <parpeoptimization/multiStartOptimization.h>
#include <parpeoptimization/optimizationProblem.h>
#include <parpeamici/amiciSimulationRunner.h>
#include <parpeoptimization/minibatchOptimization.h>

#include <amici/amici.h>
#include <amici/serialization.h>

#include <boost/serialization/map.hpp>

#include <gsl/gsl-lite.hpp>

#include <memory>
#include <cstdlib>

/** @file Interfaces between AMICI model and parPE optimization problem */

namespace parpe {

class LoadBalancerMaster;
class OptimizationResultWriter;
class MultiConditionDataProviderHDF5;
class MultiConditionDataProvider;
/**
 * @brief Run AMICI simulation for the given condition, save and return results
 * @param solver
 * @param model Model for simulation. Sensitivity
 * @param conditionIdx
 * @param jobId
 * @param dataProvider
 * @param resultWriter
 * @param logLineSearch
 * @param logger
 * @return Simulation results
 */

AmiciSimulationRunner::AmiciResultPackageSimple runAndLogSimulation(amici::Solver const &solver,
        amici::Model &model,
        int conditionIdx,
        int jobId,
        const MultiConditionDataProvider *dataProvider,
        OptimizationResultWriter *resultWriter,
        bool logLineSearch,
        Logger *logger);

/**
 * @brief Run simulations (no gradient) with given parameters and collect
 * model outputs
 * @param dataProvider
 * @param loadBalancer
 * @param maxSimulationsPerPackage
 * @param resultWriter
 * @param logLineSearch
 * @param parameters Model parameters for simulation
 * @param modelOutput in: some vector reference, will be resized.
 * output: Vector of double vectors containing AMICI ReturnData::y
 * (nt x ny, column-major)
 * @param logger
 * @param cpuTime
 * @return Simulation status
 */
FunctionEvaluationStatus getModelOutputs(
        MultiConditionDataProvider *dataProvider,
        LoadBalancerMaster *loadBalancer,
        int maxSimulationsPerPackage,
        OptimizationResultWriter *resultWriter,
        bool logLineSearch,
        gsl::span<const double> parameters,
        std::vector<std::vector<double> > &modelOutput,
        Logger *logger, double *cpuTime);

/**
 * @brief Callback function for LoadBalancer
 * @param dataProvider
 * @param resultWriter
 * @param logLineSearch
 * @param buffer In/out: message buffer
 * @param jobId: In: Identifier of the job (unique up to INT_MAX)
 */
void messageHandler(MultiConditionDataProvider *dataProvider,
                    OptimizationResultWriter *resultWriter,
                    bool logLineSearch,
                    std::vector<char> &buffer, int jobId);

/**
 * @brief The AmiciSummedGradientFunction class represents a cost function
 * based on simulations of an AMICI model for different datasets
 */
class AmiciSummedGradientFunction : public SummedGradientFunction<int> {

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
            OptimizationResultWriter *resultWriter);

    virtual ~AmiciSummedGradientFunction() = default;

    virtual FunctionEvaluationStatus evaluate(
            gsl::span<const double> parameters,
            int dataset,
            double &fval,
            gsl::span<double> gradient,
            Logger *logger,
            double *cpuTime) const override;

    virtual FunctionEvaluationStatus evaluate(
            gsl::span<const double> parameters,
            std::vector<int> datasets,
            double &fval,
            gsl::span<double> gradient,
            Logger *logger,
            double *cpuTime) const override;

    /**
     * @brief Number of optimization parameters
     * @return
     */
    virtual int numParameters() const override;

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
            double *cpuTime) const;

    virtual std::vector<std::vector<double>> getAllSigmas() const;

    virtual std::vector<std::vector<double>> getAllMeasurements() const;

    /**
     * @brief Callback function for LoadBalancer
     * @param buffer In/out: message buffer
     * @param jobId: In: Identifier of the job (unique up to INT_MAX)
     */
    virtual void messageHandler(std::vector<char> &buffer, int jobId) const;

    virtual amici::ParameterScaling getParameterScaling(int parameterIndex) const;

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
                               const std::vector<int> &dataIndices,
                               Logger *logger,
                               double *cpuTime) const;

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
                            double &simulationTimeInS,
                            gsl::span<const double> optimizationParameters
                            ) const;


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
            gsl::span<double> objectiveFunctionGradient,
            gsl::span<const double> parameters) const;

    void setSensitivityOptions(bool sensiRequired) const;

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
    OptimizationResultWriter *getResultWriter();

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


    int getNumberOfStarts() const override;

    bool restartOnFailure() const override;

    std::unique_ptr<OptimizationProblem> getLocalProblem(
            int multiStartIndex) const override;

private:
    MultiConditionDataProviderHDF5 *dp = nullptr;
    OptimizationOptions options;
    OptimizationResultWriter *resultWriter = nullptr;
    LoadBalancerMaster *loadBalancer = nullptr;
    std::unique_ptr<Logger> logger;
};


void saveSimulation(
        H5::H5File const& file, const std::string &pathStr,
        const std::vector<double> &parameters, double llh,
        gsl::span<const double> gradient, double timeElapsedInSeconds,
        gsl::span<const double> states,
        gsl::span<const double> stateSensi, gsl::span<const double> outputs,
        int jobId, int status, const std::string &label);


} // namespace parpe

#endif
