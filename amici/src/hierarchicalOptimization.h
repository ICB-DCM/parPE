#ifndef HIERACHICALOPTIMIZATION_H
#define HIERACHICALOPTIMIZATION_H

#include <optimizationProblem.h>
#include <multiConditionProblem.h>
#include <misc.h>
#include <parpeException.h>

#include <memory>
#include <cmath>
#include <numeric>

#include <H5Cpp.h>

namespace parpe {

// Currently using enum from amici enum class ParameterTransformation { none, log10 };
enum class ErrorModel { normal }; // TODO logNormal, laplace

class AnalyticalParameterProvider;
class AnalyticalParameterHdf5Reader;

/**
 * @brief The HierachicalOptimizationWrapper class is a wrapper for hierarchical optimization of
 * scaling parameters.
 *
 * Parameters with the given indices are hidden by the wrapper and computed analytically internally.
 */
class HierachicalOptimizationWrapper : public GradientFunction
{
public:
    HierachicalOptimizationWrapper(std::unique_ptr<AmiciSummedGradientFunction<int>> fun,
                                   int numConditions = 0,
                                   int numObservables = 0,
                                   int numTimepoints = 0);


    HierachicalOptimizationWrapper(std::unique_ptr<AmiciSummedGradientFunction<int>> fun,
                                   H5::H5File file, std::string hdf5RootPath,
                                   int numConditions,
                                   int numObservables,
                                   int numTimepoints,
                                   ErrorModel errorModel);

    HierachicalOptimizationWrapper(std::unique_ptr<AmiciSummedGradientFunction<int>> fun,
                                   std::unique_ptr<AnalyticalParameterProvider> scalingReader,
                                   std::unique_ptr<AnalyticalParameterProvider> offsetReader,
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
    virtual FunctionEvaluationStatus evaluate(
            const double* const parameters,
            double &fval,
            double* gradient) const override;

    std::vector<double> getDefaultScalingFactors() const;

    std::vector<double> getDefaultOffsetParameters() const;

    /**
     * @brief Run simulations with scaling parameters set to 1.0 and collect model outputs
     * @param reducedParameters parameter vector for `fun` without scaling parameters
     * @return Vector of double vectors containing AMICI ReturnData::y (nt x ny, column-major)
     */
    std::vector <std::vector<double>> getUnscaledModelOutputs(double const * const reducedParameters) const;


    /**
     * @brief Compute proportionality factors
     * @param modelOutputs Model outputs as provided by getModelOutputs
     * @return the computed scaling factors
     */
    std::vector<double> computeAnalyticalScalings(const std::vector<std::vector<double> > &measurements, std::vector <std::vector<double>>& modelOutputs) const;

    void applyOptimalScalings(std::vector<double> const& proportionalityFactors, std::vector<std::vector<double> > &modelOutputs) const;


    /**
     * @brief Compute offset parameters
     * @param modelOutputs Model outputs as provided by getModelOutputs
     * @return the computed offset parameters
     */
    std::vector<double> computeAnalyticalOffsets(const std::vector<std::vector<double> > &measurements, std::vector <std::vector<double>>& modelOutputs) const;

    void applyOptimalOffsets(std::vector<double> const& proportionalityFactors, std::vector<std::vector<double> > &modelOutputs) const;


    /**
     * @brief Compute the proportionality factor for the given observable.
     *
     * See Supplement 1.1 of [1].
     *
     * [1] Loos, Krause, Hasenauer. Hierarchical optimization for the efficient parametrization of ODE models.
     * @param
     * @return
     */

    virtual double computeAnalyticalScalings(int scalingIdx, std::vector <std::vector<double>> const& modelOutputsUnscaled,
                                  std::vector <std::vector<double>> const& measurements) const;

    virtual double computeAnalyticalOffsets(int offsetIdx, std::vector <std::vector<double>> const& modelOutputsUnscaled,
                                  std::vector <std::vector<double>> const& measurements) const;

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

    virtual FunctionEvaluationStatus evaluateWithOptimalParameters(const double * const reducedParameters, const std::vector<double> &scalings,
            const std::vector<double> &offsets,
            std::vector<std::vector<double> > &modelOutputsUnscaled,
            double &fval,
            double *gradient) const;

    /**
     * @brief Compute loglikelihood for normal distribution based on the model outputs and measurements for multiple conditions.
     * @param modelOutputsScaled
     * @return
     */
    double computeNegLogLikelihood(std::vector <std::vector<double>> const& measurements, std::vector <std::vector<double>> const& modelOutputsScaled) const;

    /**
     * @brief Compute loglikelihood for normal distribution based on the model outputs and measurements for a single condition.
     * @param modelOutputsScaled
     * @return
     */
    double computeNegLogLikelihood(std::vector<double> const& measurements, std::vector<double> const& modelOutputsScaled) const;

    virtual int numParameters() const override;

    int numScalingFactors() const;

    std::vector<int> const& getProportionalityFactorIndices() const;

    int numOffsetParameters() const;

    std::vector<int> const& getOffsetParameterIndices() const;

    std::vector<int> getAnalyticalParameterIndices() const;

    std::unique_ptr<AmiciSummedGradientFunction<int>> fun;

private:
    void init();

    /** Reads scaling parameter information from HDF5 file */
    std::unique_ptr<AnalyticalParameterProvider> scalingReader;
    /** Reads offset parameter information from HDF5 file */
    std::unique_ptr<AnalyticalParameterProvider> offsetReader;

    /** Sorted list of the indices of the scaling parameters
      * (sorting makes it easier to splice scaling and remaining parameters in getFullParameters) */
    std::vector<int> proportionalityFactorIndices;
    std::vector<int> offsetParameterIndices;

    /** Total number of conditions used in `fun` */
    int numConditions;
    /** Total number of observables occuring in `fun` */
    int numObservables;
    /** Total number of timepoints used in `fun` */
    int numTimepoints;

    /** Error model to use for computing analytical parameters and likelihood */
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
    AnalyticalParameterHdf5Reader(H5::H5File file,
                                  std::string analyticalParameterIndicesPath,
                                  std::string mapPath);

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
    int getNumAnalyticalParameters(H5::DataSet &dataset) const;

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
    HierachicalOptimizationProblemWrapper(std::unique_ptr<OptimizationProblem> problemToWrap, MultiConditionDataProvider const* dataProvider);

    virtual ~HierachicalOptimizationProblemWrapper();

    virtual void fillInitialParameters(double *buffer) const override;

    virtual void fillParametersMin(double *buffer) const override;

    virtual void fillParametersMax(double *buffer) const override;

    void fillFilteredParams(std::vector<double> const& fullParams, double *buffer) const;

    OptimizationOptions const& getOptimizationOptions() const override { return wrappedProblem->getOptimizationOptions(); }
    void setOptimizationOptions(OptimizationOptions const& options) override { wrappedProblem->setOptimizationOptions(options); }

    // TODO: need to ensure that this will work with the reduced number of parameters
    virtual std::unique_ptr<OptimizationReporter> getReporter() const override { return wrappedProblem->getReporter(); }

private:
    std::unique_ptr<OptimizationProblem> wrappedProblem;
};


/**
 * @brief Filter a vector using a list of exclude-indices. Write filter list to buffer.
 * @param valuesToFilter Original list
 * @param sortedIndicesToExclude Blacklist of indices
 * @param result Buffer to write the filtered list to. Must be at least of length valuesToFilter.size()-sortedIndicesToExclude.size().
 */
void fillFilteredParams(std::vector<double> const& valuesToFilter, const std::vector<int> &sortedIndicesToExclude, double *result);

double getDefaultScalingFactor(amici::AMICI_parameter_scaling scaling);

double getDefaultOffsetParameter(amici::AMICI_parameter_scaling scaling);

double computeAnalyticalScalings(int scalingIdx, amici::AMICI_parameter_scaling scale,
                                 const std::vector<std::vector<double> > &modelOutputsUnscaled,
                                 const std::vector<std::vector<double> > &measurements,
                                 AnalyticalParameterProvider& scalingReader,
                                 int numObservables, int numTimepoints);

double computeAnalyticalOffsets(int offsetIdx,
                                amici::AMICI_parameter_scaling scale,
                                const std::vector<std::vector<double> > &modelOutputsUnscaled,
                                const std::vector<std::vector<double> > &measurements,
                                AnalyticalParameterProvider& offsetReader,
                                int numObservables, int numTimepoints);

void applyOptimalScaling(int scalingIdx, double scalingLin,
                         std::vector<std::vector<double> > &modelOutputs,
                         AnalyticalParameterProvider& scalingReader,
                         int numObservables, int numTimepoints);

double getScaledParameter(double parameter, amici::AMICI_parameter_scaling scale);

double getUnscaledParameter(double parameter, amici::AMICI_parameter_scaling scale);

void applyOptimalOffset(int offsetIdx, double offsetLin,
                        std::vector<std::vector<double> > &modelOutputs,
                        AnalyticalParameterProvider& offsetReader,
                        int numObservables, int numTimepoints);

/**
 * @brief Assemble full parameter vector of wrapped problem from scaling parameters and numerically optimized parameters
 * @param reducedParameters
 * @param scalingFactors
 * @return Full parameter vector for `fun`
 */
std::vector<double> spliceParameters(const double * const reducedParameters, int numReduced,
                                     const std::vector<int> &proportionalityFactorIndices,
                                     const std::vector<int> &offsetParameterIndices,
                                     const std::vector<double> &scalingFactors,
                                     const std::vector<double> &offsetParameters);
} //namespace parpe

#endif // HIERACHICALOPTIMIZATION_H
