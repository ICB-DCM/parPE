#ifndef HIERARCHICALOPTIMIZATION_H
#define HIERARCHICALOPTIMIZATION_H

#include <parpeamici/hierarchicalOptimizationAnalyticalParameterProvider.h>
#include <parpeamici/multiConditionProblem.h>
#include <parpeoptimization/optimizationProblem.h>

#include <gsl/gsl-lite.hpp>

#include <memory>


namespace parpe {

enum class ErrorModel
{
    normal
}; // TODO logNormal, Laplace

class AnalyticalParameterProvider;
class AnalyticalParameterHdf5Reader;
class HierarchicalOptimizationProblemWrapper;
class HierarchicalOptimizationWrapper;

/**
 * @brief The HierarchicalOptimizationWrapper class is a wrapper for
 * hierarchical optimization of scaling parameters.
 *
 * Parameters with the given indices are hidden by the wrapper and computed
 * analytically internally.
 *
 * Computes the negative log likelihood for normally distributed measurement
 * (others to be added).
 */
class HierarchicalOptimizationWrapper : public GradientFunction
{
  public:
    /**
     * @brief For testing
     * @param fun
     * @param numConditions
     * @param numObservables
     * @param numTimepoints
     */
    HierarchicalOptimizationWrapper(
        AmiciSummedGradientFunction *wrapped_function,
        int numConditions = 0,
        int numObservables = 0);

    /**
     * @brief Get information on analytically computed parameters from HDF5 file
     * @param fun
     * @param file
     * @param hdf5RootPath
     * @param numConditions
     * @param numObservables
     * @param errorModel
     */
    HierarchicalOptimizationWrapper(
        AmiciSummedGradientFunction *wrapped_function,
        const H5::H5File& file,
        const std::string& hdf5RootPath,
        int numConditions,
        int numObservables,
        ErrorModel errorModel);

    /**
     * @brief Get information on analytically computed parameters from the
     * provided objects.
     * @param fun
     * @param scalingReader
     * @param offsetReader
     * @param numConditions
     * @param numObservables
     * @param errorModel
     */
    HierarchicalOptimizationWrapper(
        AmiciSummedGradientFunction *wrapped_function,
        std::unique_ptr<AnalyticalParameterProvider> scalingReader,
        std::unique_ptr<AnalyticalParameterProvider> offsetReader,
        std::unique_ptr<AnalyticalParameterProvider> sigmaReader,
        int numConditions,
        int numObservables,
        ErrorModel errorModel);

    using GradientFunction::evaluate;

    FunctionEvaluationStatus evaluate(gsl::span<double const> parameters,
                                      double& fval,
                                      gsl::span<double> gradient,
                                      Logger* logger,
                                      double* cpuTime) const override;

    FunctionEvaluationStatus evaluate(gsl::span<double const> reducedParameters,
                                      double& fval,
                                      gsl::span<double> gradient,
                                      std::vector<double>& fullParameters,
                                      std::vector<double>& fullGradient,
                                      Logger* logger,
                                      double* cpuTime) const;

    /**
     * @brief Get parameters for initial function evaluation
     * @return
     */
    [[nodiscard]] std::vector<double> getDefaultScalingFactors() const;

    /**
     * @brief Get parameters for initial function evaluation
     * @return
     */
    [[nodiscard]] std::vector<double> getDefaultOffsetParameters() const;

    [[nodiscard]] std::vector<double> getDefaultSigmaParameters() const;

    /**
     * @brief Run simulations with scaling parameters set to 1.0 and collect
     * model outputs
     * @param reducedParameters parameter vector for `fun` without scaling
     * parameters
     * @return Vector of double vectors containing AMICI ReturnData::y (nt x ny,
     * column-major)
     */
    [[nodiscard]] std::tuple<std::vector<std::vector<double>>,
               std::vector<std::vector<double>>>
    getUnscaledModelOutputsAndSigmas(
        const gsl::span<double const> reducedParameters,
        Logger* logger,
        double* cpuTime) const;

    /**
     * @brief Compute proportionality factors
     * @param modelOutputs Model outputs as provided by getModelOutputs
     * @return the computed scaling factors
     */
    [[nodiscard]] std::vector<double> computeAnalyticalScalings(
        std::vector<std::vector<double>> const& measurements,
        std::vector<std::vector<double>> const& modelOutputsUnscaled) const;

    void applyOptimalScalings(
        std::vector<double> const& proportionalityFactors,
        std::vector<std::vector<double>>& modelOutputs) const;

    /**
     * @brief Compute offset parameters
     * @param modelOutputs Model outputs as provided by getModelOutputs
     * @return the computed offset parameters
     */
    [[nodiscard]] std::vector<double> computeAnalyticalOffsets(
        const std::vector<std::vector<double>>& measurements,
        std::vector<std::vector<double>>& modelOutputsUnscaled) const;

    [[nodiscard]] std::vector<double> computeAnalyticalSigmas(
        std::vector<std::vector<double>> const& measurements,
        const std::vector<std::vector<double>>& modelOutputsScaled) const;

    void applyOptimalOffsets(
        std::vector<double> const& offsetParameters,
        std::vector<std::vector<double>>& modelOutputs) const;

    /**
     * @brief Create vector with sigma matrix for each condition and timepoints
     * from the given analytically computed sigmas
     * @param sigmas
     * @return
     */
    void fillInAnalyticalSigmas(
        std::vector<std::vector<double>>& allSigmas,
        const std::vector<double>& analyticalSigmas) const;

    /**
     * @brief Evaluate `fun` using the computed optimal scaling and offset
     * parameters.
     * @param reducedParameters Parameter vector without scaling and offset
     * parameters
     * @param scalings Optimal scaling parameters
     * @param offsets Optimal offset parameters
     * @param modelOutputsUnscaled Model outputs before applying optimal offset
     * and scaling parameters
     * @param fval out: computed function value
     * @param gradient out: computed function gradient
     * @return
     */

    virtual FunctionEvaluationStatus evaluateWithOptimalParameters(
        std::vector<double> const& fullParameters,
        std::vector<double> const& sigmas,
        std::vector<std::vector<double>> const& measurements,
        std::vector<std::vector<double>> const& modelOutputsScaled,
        std::vector<std::vector<double> > &fullSigmaMatrices,
        double& fval,
        const gsl::span<double> gradient,
        std::vector<double>& fullGradient,
        Logger* logger,
        double* cpuTime) const;

    /**
     * @brief Get number of parameters the function expects
     * @return That
     */
    [[nodiscard]] int numParameters() const override;

    [[nodiscard]] int numProportionalityFactors() const;

    [[nodiscard]] std::vector<int> const& getProportionalityFactorIndices() const;

    [[nodiscard]] int numOffsetParameters() const;

    [[nodiscard]] int numSigmaParameters() const;

    [[nodiscard]] std::vector<int> const& getOffsetParameterIndices() const;

    [[nodiscard]] std::vector<int> const& getSigmaParameterIndices() const;

    [[nodiscard]] std::vector<int> getAnalyticalParameterIndices() const;

    [[nodiscard]] AmiciSummedGradientFunction* getWrappedFunction() const;

    [[nodiscard]] std::vector<std::string> getParameterIds() const override;

  private:
    void init();
    /** Objective function of inner optimization problem */
    AmiciSummedGradientFunction *wrapped_function_;

    /** Reads scaling parameter information from HDF5 file */
    std::unique_ptr<AnalyticalParameterProvider> scalingReader;
    /** Reads offset parameter information from HDF5 file */
    std::unique_ptr<AnalyticalParameterProvider> offsetReader;
    std::unique_ptr<AnalyticalParameterProvider> sigmaReader;

    /** Sorted list of the indices of the scaling parameters
     * (sorting makes it easier to splice scaling and remaining parameters in
     * getFullParameters) */
    std::vector<int> proportionalityFactorIndices;
    std::vector<int> offsetParameterIndices;
    std::vector<int> sigmaParameterIndices;

    /** Total number of conditions used in `fun` */
    int numConditions;
    /** Total number of observables occurring in `fun` */
    int numObservables;

    /** Error model to use for computing analytical parameters and negative
     * log-likelihood */
    ErrorModel errorModel = ErrorModel::normal;
};


/**
 * @brief The HierarchicalOptimizationProblemWrapper class wraps an
 * OptimizationProblem and hides the analytically optimized parameters (from
 * starting point, parameter bounds, ...)
 *
 */
class HierarchicalOptimizationProblemWrapper : public OptimizationProblem
{
  public:
    HierarchicalOptimizationProblemWrapper() = default;

    HierarchicalOptimizationProblemWrapper(
        std::unique_ptr<OptimizationProblem> problemToWrap,
        const MultiConditionDataProviderHDF5* dataProvider);

    HierarchicalOptimizationProblemWrapper(
        std::unique_ptr<OptimizationProblem> problemToWrap,
        std::unique_ptr<HierarchicalOptimizationWrapper> costFun,
        std::unique_ptr<Logger> logger);

    HierarchicalOptimizationProblemWrapper(
        HierarchicalOptimizationProblemWrapper const& other) = delete;

    void fillInitialParameters(gsl::span<double> buffer) const override;

    void fillParametersMin(gsl::span<double> buffer) const override;

    void fillParametersMax(gsl::span<double> buffer) const override;

    void fillFilteredParams(std::vector<double> const& fullParams,
                            gsl::span<double> buffer) const;

    OptimizationOptions const& getOptimizationOptions() const override
    {
        return wrapped_problem_->getOptimizationOptions();
    }
    void setOptimizationOptions(OptimizationOptions const& options) override
    {
        wrapped_problem_->setOptimizationOptions(options);
    }

    // TODO: need to ensure that this will work with the reduced number of
    // parameters
    std::unique_ptr<OptimizationReporter> getReporter() const override;

  private:
    std::unique_ptr<OptimizationProblem> wrapped_problem_;
};

/**
 * @brief The HierarchicalOptimizationReporter class saves optimization
 * parameters of the inner optimization problem on each function evaluation
 * which would be hidden from the (outer) optimizer otherwise.
 */
class HierarchicalOptimizationReporter : public OptimizationReporter
{
  public:
    HierarchicalOptimizationReporter(
        HierarchicalOptimizationWrapper* gradFun,
        std::unique_ptr<OptimizationResultWriter> rw,
        std::unique_ptr<Logger> logger);

    using GradientFunction::evaluate;

    FunctionEvaluationStatus evaluate(gsl::span<const double> parameters,
                                      double& fval,
                                      gsl::span<double> gradient,
                                      Logger* logger = nullptr,
                                      double* cpuTime = nullptr) const override;

    // bool starting(gsl::span<const double> initialParameters) const override;

    // TODO: always update final parameters
    bool iterationFinished(
        gsl::span<const double> parameters,
        double objectiveFunctionValue,
        gsl::span<const double> objectiveFunctionGradient) const override;

    bool afterCostFunctionCall(
        gsl::span<const double> parameters,
        double objectiveFunctionValue,
        gsl::span<double const> objectiveFunctionGradient) const override;

    void finished(double optimalCost,
                  gsl::span<const double> parameters,
                  int exitStatus) const override;

    std::vector<double> const& getFinalParameters() const override;

    HierarchicalOptimizationWrapper* hierarchical_wrapper_ = nullptr;

    /* In addition to the vectors for the outer optimization problem,
     * we also keep the complete ones.
     */
    mutable std::vector<double> cached_full_parameters_;
    mutable std::vector<double> cached_full_gradient_;
    // TODO should override other functions as well
    // TODO: in all functions, we need to check of the provided parameters or
    // function values match To cached ones, if we want to provide all together
    // to downstream
};

/**
 * @brief Filter a vector using a list of exclude-indices. Write filter list to
 * buffer.
 * @param valuesToFilter Original list
 * @param sortedIndicesToExclude Blacklist of indices
 * @param result Buffer to write the filtered list to. Must be at least of
 * length valuesToFilter.size()-sortedIndicesToExclude.size().
 */
void
fillFilteredParams(std::vector<double> const& valuesToFilter,
                   const std::vector<int>& sortedIndicesToExclude,
                   gsl::span<double> result);

/**
 * @brief Get value to use for scaling parameter during simulation prior to
 * computing optimal values.
 * @param scaling Expected scale of the parameter
 * @return default value
 */
double
getDefaultScalingFactor(amici::ParameterScaling scaling);

/**
 * @brief Get value to use for offset parameter during simulation prior to
 * computing optimal values.
 * @param scaling Expected scale of the parameter
 * @return default value
 */
double
getDefaultOffsetParameter(amici::ParameterScaling scaling);

/**
 * @brief Compute the proportionality factor for the given observable.
 *
 * See Supplement 1.1 of [1].
 *
 * [1] Loos, Krause, Hasenauer. Hierarchical optimization for the efficient
 * parametrization of ODE models.
 *
 * @param scalingIdx
 * @param modelOutputsUnscaled Unscaled model outputs
 * @param measurements Measurements
 * @param scalingReader
 * @param numObservables Number of observables
 * @return
 */
double
computeAnalyticalScalings(
    int scalingIdx,
    const std::vector<std::vector<double>>& modelOutputsUnscaled,
    const std::vector<std::vector<double>>& measurements,
    const AnalyticalParameterProvider& scalingReader,
    int numObservables);

double
computeAnalyticalOffsets(
    int offsetIdx,
    const std::vector<std::vector<double>>& modelOutputsUnscaled,
    const std::vector<std::vector<double>>& measurements,
    AnalyticalParameterProvider& offsetReader,
    int numObservables);

double
computeAnalyticalSigmas(
    int sigmaIdx,
    const std::vector<std::vector<double>>& modelOutputsScaled,
    const std::vector<std::vector<double>>& measurements,
    const AnalyticalParameterProvider& sigmaReader,
    int numObservables,
    double epsilonAbs = 1e-12,
    double epsilonRel = 0.01);

void
applyOptimalScaling(int scalingIdx,
                    double scalingLin,
                    std::vector<std::vector<double>>& modelOutputs,
                    AnalyticalParameterProvider const& scalingReader,
                    int numObservables);

void
applyOptimalOffset(int offsetIdx,
                   double offsetLin,
                   std::vector<std::vector<double>>& modelOutputs,
                   AnalyticalParameterProvider const& offsetReader,
                   int numObservables);

/**
 * @brief Assemble full parameter vector of wrapped problem from scaling
 * parameters and numerically optimized parameters
 * @param reducedParameters
 * @param scalingFactors
 * @return Full parameter vector for `fun`
 */
std::vector<double>
spliceParameters(const gsl::span<double const> reducedParameters,
                 const std::vector<int>& proportionalityFactorIndices,
                 const std::vector<int>& offsetParameterIndices,
                 const std::vector<int>& sigmaParameterIndices,
                 const std::vector<double>& scalingFactors,
                 const std::vector<double>& offsetParameters,
                 const std::vector<double>& sigmaParameters);

/**
 * @brief Remove inner parameters
 * @return Outer parameters
 */
template <typename T>
std::vector<T>
removeInnerParameters(const gsl::span<T const> allParameters,
                      const std::vector<int>& proportionalityFactorIndices,
                      const std::vector<int>& offsetParameterIndices,
                      const std::vector<int>& sigmaParameterIndices)
{
    std::vector<T> outerParameters(
        allParameters.size() - proportionalityFactorIndices.size() -
        offsetParameterIndices.size() - sigmaParameterIndices.size());

    int nextOuterIdx = 0;
    for(int idxFull = 0; idxFull < static_cast<int>(allParameters.size());
         ++idxFull) {

        // skip if current parameter is scaling/offset/sigma
        if(std::find(proportionalityFactorIndices.begin(),
                      proportionalityFactorIndices.end(), idxFull)
            != std::end(proportionalityFactorIndices))
            continue;

        if(std::find(offsetParameterIndices.begin(),
                      offsetParameterIndices.end(), idxFull)
            != std::end(offsetParameterIndices))
            continue;

        if(std::find(sigmaParameterIndices.begin(),
                      sigmaParameterIndices.end(), idxFull)
            != std::end(sigmaParameterIndices))
            continue;

        // otherwise copy
        outerParameters[nextOuterIdx] = allParameters[idxFull];
        ++nextOuterIdx;
    }

    Ensures(nextOuterIdx == static_cast<int>(outerParameters.size()));
    return outerParameters;
}

/**
 * @brief From the given parameter vector, extract outer optimization
 * parameters, as defined in the file HDF5 file `parameterFile`
 * @param fullParameters
 * @param parameterFile
 * @param parameterPath
 * @return
 */
std::vector<double>
getOuterParameters(std::vector<double> const& fullParameters,
                   H5::H5File const& parameterFile,
                   std::string const& parameterPath);


/**
 * @brief Compute negative log-likelihood for normal distribution based on the
 * model outputs and measurements for multiple conditions.
 * @param measurements
 * @param modelOutputsScaled
 * @param sigmas
 * @return
 */
double
computeNegLogLikelihood(
    std::vector<std::vector<double>> const& measurements,
    std::vector<std::vector<double>> const& modelOutputsScaled,
    const std::vector<std::vector<double>>& sigmas);

/**
 * @brief Compute negative log-likelihood for normal distribution based on the
 * model outputs and measurements for a single condition.
 * @param measurements
 * @param modelOutputsScaled
 * @param sigmas
 * @return Negative log-likelihood for the given measurements and simulations,
 * assuming independently normally distributed noise
 */
double
computeNegLogLikelihood(std::vector<double> const& measurements,
                        std::vector<double> const& modelOutputsScaled,
                        const std::vector<double>& sigmas);

/**
 * @brief If sensitivities are computed w.r.t. analytically computed parameters
 * (which is unnecessary), this function checks they are below the given
 * threshold.
 * @param gradient
 * @param analyticalIndices
 * @param threshold
 */
void
checkGradientForAnalyticalParameters(std::vector<double> const& gradient,
                                     std::vector<int> const& analyticalIndices,
                                     double threshold);

} // namespace parpe

#endif // HIERARCHICALOPTIMIZATION_H
