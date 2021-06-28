#ifndef HIERARCHICALOPTIMIZATIONANALYTICALPARAMETERPROVIDER_H
#define HIERARCHICALOPTIMIZATIONANALYTICALPARAMETERPROVIDER_H

#include <vector>
#include <map>

#include <H5Cpp.h>

namespace parpe {


/**
 * @brief The AnalyticalParameterProvider class is an interface for providing
 * information on optimization parameters to be computed analytically
 * (proportionality factors, offsets, sigmas, ...).
 */
class AnalyticalParameterProvider
{
  public:
    virtual ~AnalyticalParameterProvider() = default;

    /**
     * @brief Get vector of condition indices for which the parameter with the
     * given index is used.
     * @param parameterIndex referring to the index in the analytical parameter
     * list in the HDF5 file
     * (*not* the optimization parameter index).
     * @return Vector of condition indices
     */
    virtual std::vector<int> getConditionsForParameter(
        int parameterIndex) const = 0;

    /**
     * @brief Get vector of observable indices for the specified condition for
     * which the specified parameter is used.
     * @param parameterIndex
     * @return
     */
    virtual std::vector<int> const& getObservablesForParameter(
        int parameterIndex,
        int conditionIdx) const = 0;

    /**
     * @brief Vector with indices of the of the analytically determined
     * parameters within the overall optimization parameter vector
     * @return
     */
    virtual std::vector<int> getOptimizationParameterIndices() const = 0;
};

class AnalyticalParameterProviderDefault : public AnalyticalParameterProvider
{
  public:
    AnalyticalParameterProviderDefault() = default;

    std::vector<int> getConditionsForParameter(
        int parameterIndex) const override;

    std::vector<int> const& getObservablesForParameter(
        int parameterIndex,
        int conditionIdx) const override;

    std::vector<int> getOptimizationParameterIndices() const override;

    // TODO private
    std::vector<std::vector<int>> conditionsForParameter;
    std::vector<int> optimizationParameterIndices;
    // x[scalingIdx][conditionIdx] -> std::vector of observableIndices
    std::vector<std::map<int, std::vector<int>>> mapping;
};

/**
 * @brief The AnalyticalParameterHdf5Reader class reads from an HDF5 file the
 * dependencies of experimental conditions and observables on parameters which
 * are to be computed analytically.
 *
 */
class AnalyticalParameterHdf5Reader : public AnalyticalParameterProvider
{
  public:
    AnalyticalParameterHdf5Reader() = default;

    /**
     * @brief AnalyticalParameterHdf5Reader
     * @param file
     * @param scalingParameterIndicesPath location in HDF5 file of the list of
     * indices of the analytically determined parameters within the overall
     * optimization parameters
     * @param mapPath path of to the dataset with the
     * parameter-observable-condition mapping
     */
    AnalyticalParameterHdf5Reader(const H5::H5File& file,
                                  std::string analyticalParameterIndicesPath,
                                  std::string mapPath);

    AnalyticalParameterHdf5Reader(AnalyticalParameterHdf5Reader const&) =
        delete;

    /**
     * @brief Get vector of condition indices for which the parameter with the
     * given index is used.
     * @param parameterIndex referring to the index in the analytical parameter
     * list in the HDF5 file
     * (*not* the optimization parameter index).
     * @return Vector of condition indices
     */
    std::vector<int> getConditionsForParameter(
        int parameterIndex) const override;

    /**
     * @brief Get vector of observable indices for the specified condition for
     * which the specified parameter is used.
     * @param parameterIndex
     * @return
     */
    std::vector<int> const& getObservablesForParameter(
        int parameterIndex,
        int conditionIdx) const override;

    /**
     * @brief Vector with indices of the of the analytically determined
     * parameters within the overall optimization parameter vector
     * @return
     */
    std::vector<int> getOptimizationParameterIndices() const override;

    ~AnalyticalParameterHdf5Reader() override;

  private:
    /**
     * @brief Get number of analytically computed parameters
     * @param dataset Read information from this dataset.
     * @return
     */
    int getNumAnalyticalParameters() const;

    /**
     * @brief Read mapping of parameter to model outputs.
     *
     * Data is expected to come as matrix/table with the following columns:
     * - scalingParameterIndex: 0-based index referring to entries in the
     *     parameter index list
     * - conditionIdx: condition index
     * - observableIdx: index of model output
     */
    void readParameterConditionObservableMappingFromFile();
    std::vector<int> readRawMap(const H5::DataSet &dataset,
                                hsize_t& nRows,
                                hsize_t& nCols) const;

    H5::H5File file;
    std::string rootPath;
    std::string mapPath;
    std::string analyticalParameterIndicesPath;

    // x[scalingIdx][conditionIdx] -> std::vector of observableIndices
    std::vector<std::map<int, std::vector<int>>> mapping;
};


}
#endif // HIERARCHICALOPTIMIZATIONANALYTICALPARAMETERPROVIDER_H
