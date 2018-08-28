#ifndef HIERACHICALOPTIMIZATION_H
#define HIERACHICALOPTIMIZATION_H

#include <optimizationProblem.h>
#include "multiConditionProblem.h"
#include <misc.h>
#include <parpeException.h>

#include <gsl/gsl-lite.hpp>

#include <memory>
#include <cmath>
#include <numeric>

#include <H5Cpp.h>

namespace parpe {

// Currently using enum from amici enum class ParameterTransformation { none, log10 };

enum class ErrorModel { normal }; // TODO logNormal, laplace

class AnalyticalParameterProvider;
class AnalyticalParameterHdf5Reader;
class HierachicalOptimizationProblemWrapper;
class HierachicalOptimizationWrapper;

/**
 * @brief The HierachicalOptimizationWrapper class is a wrapper for hierarchical optimization of
 * scaling parameters.
 *
 * Parameters with the given indices are hidden by the wrapper and computed analytically internally.
 *
 * Computes the negative log likelihood for normally distributed measurement (others to be added).
 */
class HierachicalOptimizationWrapper : public GradientFunction
{
public:
    /**
     * @brief For testing
     * @param fun
     * @param numConditions
     * @param numObservables
     * @param numTimepoints
     */
    HierachicalOptimizationWrapper(std::unique_ptr<AmiciSummedGradientFunction<int>> fun,
                                   int numConditions = 0,
                                   int numObservables = 0,
                                   int numTimepoints = 0);

    /**
     * @brief Get information on analytically computed parameters from HDF5 file
     * @param fun
     * @param file
     * @param hdf5RootPath
     * @param numConditions
     * @param numObservables
     * @param numTimepoints
     * @param errorModel
     */
    HierachicalOptimizationWrapper(std::unique_ptr<AmiciSummedGradientFunction<int>> fun,
                                   const H5::H5File &file,
                                   const std::string &hdf5RootPath,
                                   int numConditions,
                                   int numObservables,
                                   int numTimepoints,
                                   ErrorModel errorModel);

    /**
     * @brief Get information on analytically computed parameters from the provided objects.
     * @param fun
     * @param scalingReader
     * @param offsetReader
     * @param numConditions
     * @param numObservables
     * @param numTimepoints
     * @param errorModel
     */
    HierachicalOptimizationWrapper(std::unique_ptr<AmiciSummedGradientFunction<int>> fun,
                                   std::unique_ptr<AnalyticalParameterProvider> scalingReader,
                                   std::unique_ptr<AnalyticalParameterProvider> offsetReader,
                                   std::unique_ptr<AnalyticalParameterProvider> sigmaReader,
                                   int numConditions,
                                   int numObservables,
                                   int numTimepoints,
                                   ErrorModel errorModel);

    /**
     * @brief See base class
     * @param parameters
     * @param fval
     * @param gradient
     * @return
     */
    FunctionEvaluationStatus evaluate(
            gsl::span<double const> parameters,
            double &fval,
            gsl::span<double> gradient,
            Logger *logger,
            double *cpuTime) const override;

    FunctionEvaluationStatus evaluate(
            gsl::span<double const> reducedParameters,
            double &fval,
            gsl::span<double> gradient,
            std::vector<double> &fullParameters,
            std::vector<double> &fullGradient,
            Logger *logger,
            double *cpuTime) const;

    /**
     * @brief Get parameters for initial function evaluation
     * @return
     */
    std::vector<double> getDefaultScalingFactors() const;

    /**
     * @brief Get parameters for initial function evaluation
     * @return
     */
    std::vector<double> getDefaultOffsetParameters() const;

    std::vector<double> getDefaultSigmaParameters() const;

    /**
     * @brief Run simulations with scaling parameters set to 1.0 and collect model outputs
     * @param reducedParameters parameter vector for `fun` without scaling parameters
     * @return Vector of double vectors containing AMICI ReturnData::y (nt x ny, column-major)
     */
    std::vector <std::vector<double>> getUnscaledModelOutputs(const gsl::span<double const> reducedParameters, Logger *logger, double *cpuTime) const;


    /**
     * @brief Compute proportionality factors
     * @param modelOutputs Model outputs as provided by getModelOutputs
     * @return the computed scaling factors
     */
    std::vector<double> computeAnalyticalScalings(std::vector<std::vector<double>> const& measurements,
                                                  std::vector<std::vector<double>> const& modelOutputsUnscaled) const;

    void applyOptimalScalings(std::vector<double> const& proportionalityFactors,
                              std::vector<std::vector<double> > &modelOutputs) const;


    /**
     * @brief Compute offset parameters
     * @param modelOutputs Model outputs as provided by getModelOutputs
     * @return the computed offset parameters
     */
    std::vector<double> computeAnalyticalOffsets(const std::vector<std::vector<double> > &measurements,
                                                 std::vector <std::vector<double>>& modelOutputsUnscaled) const;

    std::vector<double> computeAnalyticalSigmas( std::vector<std::vector<double> > const& measurements,
                                                 std::vector <std::vector<double>>& modelOutputsScaled) const;

    void applyOptimalOffsets(std::vector<double> const& offsetParameters,
                             std::vector<std::vector<double> > &modelOutputs) const;

    /**
     * @brief Create vector with sigma matrix for each condition and timepoints from the given analytically computed sigmas
     * @param sigmas
     * @return
     */
    void fillInAnalyticalSigmas(std::vector<std::vector<double> > &allSigmas,
                                const std::vector<double> &analyticalSigmas) const;

    /**
     * @brief Evaluate `fun` using the computed optimal scaling and offset parameters.
     * @param reducedParameters Parameter vector without scaling and offset parameters
     * @param scalings Optimal scaling parameters
     * @param offsets Optimal offset parameters
     * @param modelOutputsUnscaled Model outputs before applying optimal offset and scaling parameters
     * @param fval out: computed function value
     * @param gradient out: computed function gradient
     * @return
     */

    virtual FunctionEvaluationStatus evaluateWithOptimalParameters(
            std::vector<double> const& fullParameters,
            std::vector<double> const& sigmas,
            std::vector<std::vector<double>> const& measurements,
            std::vector<std::vector<double>> const& modelOutputsScaled,
            double &fval,
            const gsl::span<double> gradient, std::vector<double> &fullGradient,
            Logger *logger, double *cpuTime) const;

    /**
     * @brief Get number of parameters the function expects
     * @return That
     */
    virtual int numParameters() const override;

    int numProportionalityFactors() const;

    std::vector<int> const& getProportionalityFactorIndices() const;

    int numOffsetParameters() const;

    int numSigmaParameters() const;

    std::vector<int> const& getOffsetParameterIndices() const;

    std::vector<int> const& getSigmaParameterIndices() const;

    std::vector<int> getAnalyticalParameterIndices() const;

    std::unique_ptr<AmiciSummedGradientFunction<int>> fun;

private:
    void init();

    /** Reads scaling parameter information from HDF5 file */
    std::unique_ptr<AnalyticalParameterProvider> scalingReader;
    /** Reads offset parameter information from HDF5 file */
    std::unique_ptr<AnalyticalParameterProvider> offsetReader;
    std::unique_ptr<AnalyticalParameterProvider> sigmaReader;

    /** Sorted list of the indices of the scaling parameters
      * (sorting makes it easier to splice scaling and remaining parameters in getFullParameters) */
    std::vector<int> proportionalityFactorIndices;
    std::vector<int> offsetParameterIndices;
    std::vector<int> sigmaParameterIndices;

    /** Total number of conditions used in `fun` */
    int numConditions;
    /** Total number of observables occuring in `fun` */
    int numObservables;
    /** Total number of timepoints used in `fun` */
    int numTimepoints;

    /** Error model to use for computing analytical parameters and negative log-likelihood */
    ErrorModel errorModel = ErrorModel::normal;
};


/**
 * @brief The AnalyticalParameterProvider class is an interface for providing information on
 * optimization parameters to be computed analytically (proportionality factors, offsets, sigmas, ...).
 */
class AnalyticalParameterProvider {
public:

    virtual ~AnalyticalParameterProvider() {}

    /**
     * @brief Get vector of condition indices for which the parameter with the given index is used.
     * @param parameterIndex referring to the index in the analytical parameter list in the hdf5 file
     * (*not* the optimization parameter index).
     * @return Vector of condition indice
     */
    virtual std::vector<int> getConditionsForParameter(int parameterIndex) const = 0;

    /**
     * @brief Get vector of observable indices for the specified condition for which the specified parameter is used.
     * @param parameterIndex
     * @return
     */
    virtual std::vector<int> const& getObservablesForParameter(int parameterIndex, int conditionIdx) const = 0;

    /**
     * @brief Vector with indices of the of the analytically determined parameters within the
     * overall optimization parameter vector
     * @return
     */
    virtual std::vector<int> getOptimizationParameterIndices() const = 0;

};



class AnalyticalParameterProviderDefault : public AnalyticalParameterProvider {
public:
    AnalyticalParameterProviderDefault() = default;

    std::vector<int> getConditionsForParameter(int parameterIndex) const override;

    std::vector<int> const& getObservablesForParameter(int parameterIndex, int conditionIdx) const override;

    std::vector<int> getOptimizationParameterIndices() const override;

    // TODO private
    std::vector <std::vector<int>> conditionsForParameter;
    std::vector <int> optimizationParameterIndices;
    // x[scalingIdx][conditionIdx] -> std::vector of observableIndicies
    std::vector<std::map<int, std::vector<int>>> mapping;
};


/**
 * @brief The AnalyticalParameterHdf5Reader class reads from an HDF5 file the dependencies of experimental conditions
 * and observables on parameters which are to be computed analytically.
 *
 */
class AnalyticalParameterHdf5Reader : public AnalyticalParameterProvider {
public:
    AnalyticalParameterHdf5Reader() = default;

    /**
     * @brief AnalyticalParameterHdf5Reader
     * @param file
     * @param scalingParameterIndicesPath location in hdf5 file of the list of indices
     * of the analytically determined parameters within the overall optimization parameters
     * @param mapPath path of to the dataset with the parameter-oberservable-condition mapping
     */
    AnalyticalParameterHdf5Reader(const H5::H5File &file,
                                  const std::string &analyticalParameterIndicesPath,
                                  const std::string &mapPath);

    AnalyticalParameterHdf5Reader(AnalyticalParameterHdf5Reader const&) = delete;

    /**
     * @brief Get vector of condition indices for which the parameter with the given index is used.
     * @param parameterIndex referring to the index in the analytical parameter list in the hdf5 file
     * (*not* the optimization parameter index).
     * @return Vector of condition indice
     */
    std::vector<int> getConditionsForParameter(int parameterIndex) const override;

    /**
     * @brief Get vector of observable indices for the specified condition for which the specified parameter is used.
     * @param parameterIndex
     * @return
     */
    std::vector<int> const& getObservablesForParameter(int parameterIndex, int conditionIdx) const override;

    /**
     * @brief Vector with indices of the of the analytically determined parameters within the
     * overall optimization parameter vector
     * @return
     */
    std::vector<int> getOptimizationParameterIndices() const override;

private:
    /**
     * @brief Get number of analytically computed parameters
     * @param dataset Read information from this dataset.
     * @return
     */
    int getNumAnalyticalParameters() const;

    void readParameterConditionObservableMappingFromFile();
    std::vector<int> readRawMap(H5::DataSet& dataset, hsize_t &nRows, hsize_t &nCols);

    H5::H5File file;
    std::string rootPath;
    std::string mapPath;
    std::string analyticalParameterIndicesPath;

    // x[scalingIdx][conditionIdx] -> std::vector of observableIndicies
    std::vector<std::map<int, std::vector<int>>> mapping;
};


/**
 * @brief The HierachicalOptimizationProblemWrapper class wraps an OptimizationProblem
 * and hides the analytically optimizated parameters (from starting point, parameter bounds, ...)
 *
 */
class HierachicalOptimizationProblemWrapper : public OptimizationProblem {
public:
    HierachicalOptimizationProblemWrapper() = default;

    HierachicalOptimizationProblemWrapper(std::unique_ptr<OptimizationProblem> problemToWrap,
                                          const MultiConditionDataProviderHDF5 *dataProvider);

    HierachicalOptimizationProblemWrapper(std::unique_ptr<OptimizationProblem> problemToWrap,
                                          std::unique_ptr<HierachicalOptimizationWrapper> costFun,
                                          std::unique_ptr<Logger> logger);

    HierachicalOptimizationProblemWrapper(HierachicalOptimizationProblemWrapper const& other) = delete;

    virtual ~HierachicalOptimizationProblemWrapper() override;

    virtual void fillInitialParameters(gsl::span<double> buffer) const override;

    virtual void fillParametersMin(gsl::span<double> buffer) const override;

    virtual void fillParametersMax(gsl::span<double> buffer) const override;

    void fillFilteredParams(std::vector<double> const& fullParams, gsl::span<double> buffer) const;

    OptimizationOptions const& getOptimizationOptions() const override { return wrappedProblem->getOptimizationOptions(); }
    void setOptimizationOptions(OptimizationOptions const& options) override { wrappedProblem->setOptimizationOptions(options); }

    // TODO: need to ensure that this will work with the reduced number of parameters
    virtual std::unique_ptr<OptimizationReporter> getReporter() const override;

private:
    std::unique_ptr<OptimizationProblem> wrappedProblem;
};


/**
 * @brief The HierarchicalOptimizationReporter class saves optimization parameters
 * of the inner optimization problem on each function evaluation which would be
 * hidden from the (outer) optimizer otherwise.
 */
class HierarchicalOptimizationReporter : public OptimizationReporter {
public:
    HierarchicalOptimizationReporter(HierachicalOptimizationWrapper *gradFun,
                                     std::unique_ptr<OptimizationResultWriter> rw,
                                     std::unique_ptr<Logger> logger);

    FunctionEvaluationStatus evaluate(
            gsl::span<const double> parameters, double &fval,
            gsl::span<double> gradient,
            Logger *logger = nullptr,
            double *cpuTime = nullptr) const override;

    // bool starting(gsl::span<const double> initialParameters) const override;

    // TODO: always update final parameters
    virtual bool iterationFinished(gsl::span<const double> parameters,
                                   double objectiveFunctionValue,
                                   gsl::span<const double> objectiveFunctionGradient) const override;

    // virtual bool beforeCostFunctionCall(gsl::span<const double> parameters) const override;

    virtual bool afterCostFunctionCall(gsl::span<const double> parameters,
                                       double objectiveFunctionValue,
                                       gsl::span<double const> objectiveFunctionGradient) const override;

    void finished(double optimalCost, gsl::span<const double> parameters, int exitStatus) const override;

    std::vector<double> const& getFinalParameters() const override;

    HierachicalOptimizationWrapper *hierarchicalWrapper = nullptr;

    /* In addition to the vectors for the outer optimization problem,
     * we also keep the complete ones.
     */
    mutable std::vector<double> cachedFullParameters;
    mutable std::vector<double> cachedFullGradient;
    // TODO should override other functions as well
    // TODO: in all functions, we need to check of the provided parameters or functio nvalues match
    // To cached ones, if we want to provide all together to downstreams

};

/**
 * @brief Filter a vector using a list of exclude-indices. Write filter list to buffer.
 * @param valuesToFilter Original list
 * @param sortedIndicesToExclude Blacklist of indices
 * @param result Buffer to write the filtered list to. Must be at least of length valuesToFilter.size()-sortedIndicesToExclude.size().
 */
void fillFilteredParams(std::vector<double> const& valuesToFilter,
                        const std::vector<int> &sortedIndicesToExclude,
                        gsl::span<double> result);

double getDefaultScalingFactor(amici::AMICI_parameter_scaling scaling);

double getDefaultOffsetParameter(amici::AMICI_parameter_scaling scaling);

/**
 * @brief Compute the proportionality factor for the given observable.
 *
 * See Supplement 1.1 of [1].
 *
 * [1] Loos, Krause, Hasenauer. Hierarchical optimization for the efficient parametrization of ODE models.
 * @param
 * @return
 */

double computeAnalyticalScalings(int scalingIdx,
                                 const std::vector<std::vector<double> > &modelOutputsUnscaled,
                                 const std::vector<std::vector<double> > &measurements,
                                 const AnalyticalParameterProvider &scalingReader,
                                 int numObservables, int numTimepoints);

double computeAnalyticalOffsets(int offsetIdx,
                                const std::vector<std::vector<double> > &modelOutputsUnscaled,
                                const std::vector<std::vector<double> > &measurements,
                                AnalyticalParameterProvider& offsetReader,
                                int numObservables, int numTimepoints);

double computeAnalyticalSigmas(int sigmaIdx,
                               const std::vector<std::vector<double> > &modelOutputsScaled,
                               const std::vector<std::vector<double> > &measurements,
                               const AnalyticalParameterProvider &sigmaReader,
                               int numObservables, int numTimepoints);

void applyOptimalScaling(int scalingIdx, double scalingLin,
                         std::vector<std::vector<double> > &modelOutputs,
                         AnalyticalParameterProvider const& scalingReader,
                         int numObservables, int numTimepoints);

double getScaledParameter(double parameter, amici::AMICI_parameter_scaling scale);

double getUnscaledParameter(double parameter, amici::AMICI_parameter_scaling scale);

void applyOptimalOffset(int offsetIdx, double offsetLin,
                        std::vector<std::vector<double> > &modelOutputs,
                        AnalyticalParameterProvider const& offsetReader,
                        int numObservables, int numTimepoints);

/**
 * @brief Assemble full parameter vector of wrapped problem from scaling parameters and numerically optimized parameters
 * @param reducedParameters
 * @param scalingFactors
 * @return Full parameter vector for `fun`
 */
std::vector<double> spliceParameters(const gsl::span<double const> reducedParameters,
                                     const std::vector<int> &proportionalityFactorIndices,
                                     const std::vector<int> &offsetParameterIndices,
                                     const std::vector<int> &sigmaParameterIndices,
                                     const std::vector<double> &scalingFactors,
                                     const std::vector<double> &offsetParameters,
                                     const std::vector<double> &sigmaParameters);

/**
 * @brief Compute negative log-likelihood for normal distribution based on the model outputs and measurements for multiple conditions.
 * @param modelOutputsScaled
 * @return
 */
double computeNegLogLikelihood(std::vector <std::vector<double>> const& measurements,
                               std::vector <std::vector<double>> const& modelOutputsScaled,
                               const std::vector<std::vector<double> > &sigmas);

/**
 * @brief Compute negative log-likelihood for normal distribution based on the model outputs and measurements for a single condition.
 * @param modelOutputsScaled
 * @return
 */
double computeNegLogLikelihood(std::vector<double> const& measurements,
                               std::vector<double> const& modelOutputsScaled,
                               const std::vector<double> &sigmas);

void checkGradientForAnalyticalParameters(std::vector<double> const& gradient,
                                          std::vector<int> const& analyticalIndices, double threshold);

} //namespace parpe

#endif // HIERACHICALOPTIMIZATION_H
