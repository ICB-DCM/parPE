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
                                   H5::H5File file, std::string hdf5RootPath,
                                   int numConditions,
                                   int numObservables,
                                   int numTimepoints,
                                   ErrorModel errorModel);

    HierachicalOptimizationWrapper(std::unique_ptr<AmiciSummedGradientFunction<int>> fun,
                                   std::unique_ptr<AnalyticalParameterHdf5Reader> reader,
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

    void applyOptimalScaling(int scalingIdx, double scaling, std::vector <std::vector<double>>&  modelOutputs) const;


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

    virtual FunctionEvaluationStatus evaluateWithScalings(const double* const reducedParameters, std::vector<double> const &scalings,
            std::vector<std::vector<double> > &modelOutputsUnscaled,
            double &fval,
            double* gradient) const;


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


    /**
     * @brief Assemble full parameter vector of wrapped problem from scaling parameters and numerically optimized parameters
     * @param reducedParameters
     * @param scalingFactors
     * @return Full parameter vector for `fun`
     */
    virtual std::vector<double> spliceParameters(double const * const reducedParameters, int numReduced, std::vector<double> const& scalingFactors) const;

    virtual int numParameters() const override;

    int numScalingFactors() const;

    std::vector<int> const& getProportionalityFactorIndices() const;

    std::unique_ptr<AmiciSummedGradientFunction<int>> fun;

private:
    /** Sorted list of the indices of the scaling parameters
      * (sorting makes it easier to splice scaling and remaining parameters in getFullParameters) */
    std::unique_ptr<AnalyticalParameterHdf5Reader> reader;
    std::vector<int> proportionalityFactorIndices;
    int numConditions;
    int numObservables;
    int numTimepoints;

    ErrorModel errorModel;
};


/**
 * @brief The AnalyticalParameterHdf5Reader class reads from an HDF5 file the dependencies of experimental conditions
 * and observables on parameters which are to be computed analytically.
 *
 * TODO:      * @param parameterScalingPath List stating whether the wrapped functions considers
     * these parameters as log-scaled, or ... parameters
 * This should come not from here, but from the wrapped functions above -> add method there
 */
class AnalyticalParameterHdf5Reader {
public:
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
    std::vector<int> getConditionsForParameter(int parameterIndex) const;

    /**
     * @brief Get vector of observable indices for the specified condition for which the specified parameter is used.
     * @param parameterIndex
     * @return
     */
    std::vector<int> const& getObservablesForParameter(int parameterIndex, int conditionIdx) const;

    /**
     * @brief Vector with indices of the of the analytically determined parameters within the
     * overall optimization parameter vector
     * @return
     */
    std::vector<int> getOptimizationParameterIndices();

    int getNumAnalyticalParameters(H5::DataSet &dataset) const;
private:
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
 * @brief The HierachicalOptimizationProblemWrapper class wraps an OptimizationProblem and hides the analytically optimizated parameters (from starting point, parameter bounds, ...)
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

    std::unique_ptr<OptimizationProblem> wrappedProblem;
};

void fillFilteredParams(std::vector<double> const& fullParams, const std::vector<int> &sortedFilterIndices, double *buffer);

} //namespace parpe

#endif // HIERACHICALOPTIMIZATION_H
