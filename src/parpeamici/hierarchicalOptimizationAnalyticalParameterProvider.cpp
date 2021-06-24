#include <parpeamici/hierarchicalOptimizationAnalyticalParameterProvider.h>

#include <parpecommon/hdf5Misc.h>
#include <parpecommon/parpeException.h>

namespace parpe {


std::vector<int>
AnalyticalParameterProviderDefault::getConditionsForParameter(
    int parameterIndex) const
{
    return conditionsForParameter[parameterIndex];
}

const std::vector<int>&
AnalyticalParameterProviderDefault::getObservablesForParameter(
    int parameterIndex,
    int conditionIdx) const
{
    return mapping[parameterIndex].at(conditionIdx);
}

std::vector<int>
AnalyticalParameterProviderDefault::getOptimizationParameterIndices() const
{
    return optimizationParameterIndices;
}

AnalyticalParameterHdf5Reader::AnalyticalParameterHdf5Reader(
    H5::H5File const& file,
    std::string analyticalParameterIndicesPath,
    std::string mapPath)
    : mapPath(std::move(mapPath))
      , analyticalParameterIndicesPath(std::move(analyticalParameterIndicesPath))
{
    [[maybe_unused]] auto lock = hdf5MutexGetLock();
    this->file = file; // copy while mutex is locked!
    readParameterConditionObservableMappingFromFile();
}


std::vector<int>
AnalyticalParameterHdf5Reader::getConditionsForParameter(
    int parameterIndex) const
{
    std::vector<int> result;
    result.reserve(mapping[parameterIndex].size());
    for (auto const& kvp : mapping[parameterIndex])
        result.push_back(kvp.first);
    return result;
}

const std::vector<int>&
AnalyticalParameterHdf5Reader::getObservablesForParameter(
    int parameterIndex,
    int conditionIdx) const
{
    return mapping[parameterIndex].at(conditionIdx);
}

std::vector<int>
AnalyticalParameterHdf5Reader::getOptimizationParameterIndices() const
{
    [[maybe_unused]] auto lock = hdf5MutexGetLock();
    std::vector<int> analyticalParameterIndices;

    if (file.nameExists(analyticalParameterIndicesPath)) {
        auto dataset = file.openDataSet(analyticalParameterIndicesPath);
        auto dataspace = dataset.getSpace();

        auto ndims = dataspace.getSimpleExtentNdims();
        if (ndims != 1)
            throw ParPEException(
                "Invalid dimension in getOptimizationParameterIndices.");
        hsize_t numScalings = 0;
        dataspace.getSimpleExtentDims(&numScalings);

        analyticalParameterIndices.resize(numScalings);
        dataset.read(analyticalParameterIndices.data(),
                     H5::PredType::NATIVE_INT);
    }
    return analyticalParameterIndices;
}

int
AnalyticalParameterHdf5Reader::getNumAnalyticalParameters() const
{
    hsize_t numAnalyticalParameters = 0;
    [[maybe_unused]] auto lock = hdf5MutexGetLock();

    if (file.nameExists(analyticalParameterIndicesPath)) {
        auto dataset = file.openDataSet(analyticalParameterIndicesPath);
        auto dataspace = dataset.getSpace();
        auto ndims = dataspace.getSimpleExtentNdims();
        if (ndims != 1)
            throw ParPEException(
                "Invalid dimension in getOptimizationParameterIndices.");
        dataspace.getSimpleExtentDims(&numAnalyticalParameters);
    }

    return numAnalyticalParameters;
}

void
AnalyticalParameterHdf5Reader::readParameterConditionObservableMappingFromFile()
{
    [[maybe_unused]] auto lock = hdf5MutexGetLock();
    H5_SAVE_ERROR_HANDLER;
    try {
        int numScalings = getNumAnalyticalParameters();
        auto dataset = file.openDataSet(mapPath);
        if (numScalings == 0)
            return;

        // column indices in dataspace
        constexpr int parameterCol = 0;
        constexpr int conditionCol = 1;
        constexpr int observableCol = 2;

        mapping.resize(numScalings);

        hsize_t nRows = 0, nCols = 0;
        auto rawMap = readRawMap(dataset, nRows, nCols);

        for (int i = 0; (unsigned)i < nRows; ++i) {
            int scalingIdx = rawMap[i * nCols + parameterCol];
            int conditionIdx = rawMap[i * nCols + conditionCol];
            int observableIdx = rawMap[i * nCols + observableCol];
            mapping[scalingIdx][conditionIdx].push_back(observableIdx);
        }
    } catch (H5::FileIException&) {
        return;
    }
    H5_RESTORE_ERROR_HANDLER;
}

std::vector<int>
AnalyticalParameterHdf5Reader::readRawMap(H5::DataSet const& dataset,
                                          hsize_t& nRows,
                                          hsize_t& nCols)
{
    [[maybe_unused]] auto lock = hdf5MutexGetLock();

    auto dataspace = dataset.getSpace();
    auto ndims = dataspace.getSimpleExtentNdims();
    if (ndims != 2)
        throw ParPEException(
            "Invalid dimension for analytical parameter map, expected 2.");

    hsize_t dims[ndims];
    dataspace.getSimpleExtentDims(dims);
    nRows = dims[0];
    nCols = dims[1];
    if (nRows && nCols != 3)
        throw ParPEException("Invalid dimension for analytical parameter map, "
                             "expected 3 columns.");

    std::vector<int> rawMap(nRows * nCols);
    dataset.read(rawMap.data(), H5::PredType::NATIVE_INT);

    return rawMap;
}

}
