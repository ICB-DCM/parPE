#include "hierachicalOptimization.h"

#include <assert.h>
#include <exception>

#ifdef __INTEL_COMPILER
// constexpr did not work on icc (ICC) 16.0.4 20160811
#define constexpr
#endif

namespace parpe {


HierachicalOptimizationWrapper::HierachicalOptimizationWrapper(
        std::unique_ptr<AmiciSummedGradientFunction<int> > fun,
        std::unique_ptr<parpe::ScalingFactorHdf5Reader> scalingReader,
        std::unique_ptr<parpe::OffsetParameterHdf5Reader> offsetReader,
        int numConditions, int numObservables, int numTimepoints,
        ParameterTransformation parameterTransformation,
        ErrorModel errorModel)
    : fun(std::move(fun)),
      scalingReader(std::move(scalingReader)),
      offsetReader(std::move(offsetReader)),
      numConditions(numConditions),
      numObservables(numObservables),
      numTimepoints(numTimepoints),
      parameterTransformation(parameterTransformation),
      errorModel(errorModel) {
      this->parameterTransformation = ParameterTransformation::none;
    //if(parameterTransformation != ParameterTransformation::log10 || errorModel != ErrorModel::normal) {
    //    throw ParPEException("Only log10 and gaussian noise are supported so far.");
    //}

    proportionalityFactorIndices = this->scalingReader->getProportionalityFactorIndices();
    std::sort(this->proportionalityFactorIndices.begin(),
              this->proportionalityFactorIndices.end());

    offsetParameterIndices = this->offsetReader->getOffsetParameterIndices();
    std::sort(this->offsetParameterIndices.begin(),
              this->offsetParameterIndices.end());
}

FunctionEvaluationStatus HierachicalOptimizationWrapper::evaluate(const double * const parameters, double &fval, double *gradient) const {
    if(numScalingFactors() == 0) {
        // nothing to do, just pass through
        std::vector<int> dataIndices(numConditions);
        std::iota(dataIndices.begin(), dataIndices.end(), 0);

        return fun->evaluate(parameters, dataIndices, fval, gradient);
    }

    // evaluate with scaling parameters set to 1
    auto modelOutput = getUnscaledModelOutputs(parameters);

    // compute correct scaling factors analytically
    auto scalings = computeAnalyticalScalings(fun->getAllMeasurements(), modelOutput);

    // compute correct offset parameters analytically
    auto offsets = computeAnalyticalOffsets(fun->getAllMeasurements(), modelOutput);
    // std::cout << "offsets:" << offsets;
    // evaluate with analytical scaling parameters
    auto status = evaluateWithScalings(parameters, scalings, offsets, modelOutput, fval, gradient);

    return status;
}

std::vector<double> HierachicalOptimizationWrapper::getDefaultScalingFactors() const
{
    double scaling = NAN;
    switch (parameterTransformation) {
    case ParameterTransformation::log10:
        scaling = 0.0;
        break;
    case ParameterTransformation::none:
        scaling = 1.0;
        break;
    }

    return std::vector<double>(numScalingFactors(), scaling);
}

std::vector<double> HierachicalOptimizationWrapper::getDefaultOffsetParameters() const
{
    double offset = NAN;
    switch (parameterTransformation) {
    case ParameterTransformation::log10:
        offset = -INFINITY;
        break;
    case ParameterTransformation::none:
        offset = 0.0;
        break;
    }

    return std::vector<double>(numOffsetParameters(), offset);
}

std::vector<std::vector<double> > HierachicalOptimizationWrapper::getUnscaledModelOutputs(const double * const reducedParameters) const {
    // run simulations (no gradient!) with scaling parameters == 1, collect outputs

    auto scalingDummy = getDefaultScalingFactors();
    auto offsetDummy = getDefaultOffsetParameters();
    // splice hidden scaling parameter and external parameters
    auto fullParameters = spliceParameters(reducedParameters, numParameters(), scalingDummy, offsetDummy);

    std::vector<std::vector<double> > modelOutput(numConditions);
    fun->getModelOutputs(fullParameters.data(), modelOutput);

    return modelOutput;
}

// Scaling code

std::vector<double> HierachicalOptimizationWrapper::computeAnalyticalScalings(std::vector<std::vector<double>> const& measurements, std::vector<std::vector<double> > &modelOutputs) const {
    // NOTE: does not handle replicates, assumes normal distribution, does not compute sigmas

    int numProportionalityFactors = proportionalityFactorIndices.size();
    std::vector<double> proportionalityFactors(numProportionalityFactors);

    for(int i = 0; i < numProportionalityFactors; ++i) {
        proportionalityFactors[i] = computeAnalyticalScalings(i, modelOutputs, measurements);
    }

    return proportionalityFactors;
}

void HierachicalOptimizationWrapper::applyOptimalScalings(std::vector<double> const& proportionalityFactors, std::vector<std::vector<double> > &modelOutputs) const {

    for(int i = 0; (unsigned) i < proportionalityFactors.size(); ++i) {
        double scaling = proportionalityFactors[i];

        switch(parameterTransformation) {
        case ParameterTransformation::none:
            break;
        case ParameterTransformation::log10:
            scaling = pow(10, scaling);
            break;
        }
        applyOptimalScaling(i, scaling, modelOutputs);
    }
}


void HierachicalOptimizationWrapper::applyOptimalScaling(int scalingIdx, double scaling, std::vector<std::vector<double> > &modelOutputs) const {
    auto dependentConditions = scalingReader->getConditionsForScalingParameter(scalingIdx);
    for (auto const conditionIdx: dependentConditions) {
        auto dependentObservables = scalingReader->getObservablesForScalingParameter(scalingIdx, conditionIdx);
        for(auto const observableIdx: dependentObservables) {
            assert(observableIdx < numObservables);
            for(int timeIdx = 0; timeIdx < numTimepoints; ++timeIdx) {
                modelOutputs[conditionIdx][observableIdx + timeIdx * numObservables] *= scaling;
            }
        }
    }
}

double HierachicalOptimizationWrapper::computeAnalyticalScalings(int scalingIdx, const std::vector<std::vector<double> > &modelOutputsUnscaled, const std::vector<std::vector<double> > &measurements) const {
    auto dependentConditions = scalingReader->getConditionsForScalingParameter(scalingIdx);

    double enumerator = 0.0;
    double denominator = 0.0;

    for (auto const conditionIdx: dependentConditions) {
        auto dependentObservables = scalingReader->getObservablesForScalingParameter(scalingIdx, conditionIdx);
        for(auto const observableIdx: dependentObservables) {
            for(int timeIdx = 0; timeIdx < numTimepoints; ++timeIdx) {
                double mes = measurements[conditionIdx][observableIdx + timeIdx * numObservables];
                if(!std::isnan(mes)) {
                    double sim = modelOutputsUnscaled[conditionIdx][observableIdx + timeIdx * numObservables];
                    assert(!std::isnan(sim));
                    enumerator += sim * mes;
                    denominator += sim * sim;
                }
            }
        }
    }

    double proportionalityFactor = enumerator / denominator;

    switch (parameterTransformation) {
    case ParameterTransformation::log10:
        proportionalityFactor = log10(proportionalityFactor);
        break;
    case ParameterTransformation::none:
        break;
    }

    // TODO ensure positivity!
    return proportionalityFactor;
}

// Offset code


std::vector<double> HierachicalOptimizationWrapper::computeAnalyticalOffsets(std::vector<std::vector<double>> const& measurements, std::vector<std::vector<double> > &modelOutputs) const {
    // NOTE: does not handle replicates, assumes normal distribution, does not compute sigmas

    int numOffsetParameters = offsetParameterIndices.size();
    std::vector<double> offsetParameters(numOffsetParameters);

    for(int i = 0; i < numOffsetParameters; ++i) {
        offsetParameters[i] = computeAnalyticalOffsets(i, modelOutputs, measurements);
    }

    return offsetParameters;
}

void HierachicalOptimizationWrapper::applyOptimalOffsets(std::vector<double> const& offsetParameters, std::vector<std::vector<double> > &modelOutputs) const {

    for(int i = 0; (unsigned) i < offsetParameters.size(); ++i) {
        double offset = offsetParameters[i];

        switch(parameterTransformation) {
        case ParameterTransformation::none:
            break;
        case ParameterTransformation::log10:
            offset = pow(10, offset);
            break;
        }
        applyOptimalOffset(i, offset, modelOutputs);
    }
}


void HierachicalOptimizationWrapper::applyOptimalOffset(int offsetIdx, double offset, std::vector<std::vector<double> > &modelOutputs) const {
    auto dependentConditions = offsetReader->getConditionsForOffsetParameter(offsetIdx);
    for (auto const conditionIdx: dependentConditions) {
        auto dependentObservables = offsetReader->getObservablesForOffsetParameter(offsetIdx, conditionIdx);
        for(auto const observableIdx: dependentObservables) {
            assert(observableIdx < numObservables);
            for(int timeIdx = 0; timeIdx < numTimepoints; ++timeIdx) {
                modelOutputs[conditionIdx][observableIdx + timeIdx * numObservables] += offset;
            }
        }
    }
}

double HierachicalOptimizationWrapper::computeAnalyticalOffsets(int offsetIdx, const std::vector<std::vector<double> > &modelOutputsUnscaled, const std::vector<std::vector<double> > &measurements) const {
    auto dependentConditions = offsetReader->getConditionsForOffsetParameter(offsetIdx);

    double enumerator = 0.0;
    double denominator = 0.0;

    for (auto const conditionIdx: dependentConditions) {
        auto dependentObservables = offsetReader->getObservablesForOffsetParameter(offsetIdx, conditionIdx);
        for(auto const observableIdx: dependentObservables) {
            for(int timeIdx = 0; timeIdx < numTimepoints; ++timeIdx) {
                double mes = measurements[conditionIdx][observableIdx + timeIdx * numObservables];
                if(!std::isnan(mes)) {
                    double sim = modelOutputsUnscaled[conditionIdx][observableIdx + timeIdx * numObservables];
                    assert(!std::isnan(sim));
                    enumerator += mes - sim;
                    denominator += 1;
                }
            }
        }
    }

    double offsetParameter = enumerator / denominator;

    switch (parameterTransformation) {
    case ParameterTransformation::log10:
        offsetParameter = log10(offsetParameter);
        break;
    case ParameterTransformation::none:
        break;
    }

    // TODO ensure positivity!
    return offsetParameter;
}

//

FunctionEvaluationStatus HierachicalOptimizationWrapper::evaluateWithScalings(
        // adapt to offsets
        const double * const reducedParameters,
        const std::vector<double> &scalings,
        const std::vector<double> &offsets,
        std::vector<std::vector<double> > &modelOutputsUnscaled,
        double &fval, double *gradient) const {

    if(gradient) {
        // simulate with updated theta for sensitivities
        auto fullParameters = spliceParameters(reducedParameters, numParameters(), scalings, offsets);
        // TODO: only those necessary?
        std::vector<int> dataIndices(numConditions);
        std::iota(dataIndices.begin(), dataIndices.end(), 0);

        // Need intermediary buffer because optimizer expects fewer parameters
        std::vector<double> fullGradient(fullParameters.size());
        auto status = fun->evaluate(fullParameters.data(), dataIndices, fval, fullGradient.data());
        if(status != functionEvaluationSuccess)
            return status;
        auto scalingParameterIndices = proportionalityFactorIndices;
        scalingParameterIndices.insert(scalingParameterIndices.end(),offsetParameterIndices.begin(),offsetParameterIndices.end());
        std::sort(scalingParameterIndices.begin(),scalingParameterIndices.end());
        fillFilteredParams(fullGradient, scalingParameterIndices, gradient);

    } else {
        applyOptimalScalings(scalings, modelOutputsUnscaled);
        applyOptimalOffsets(offsets, modelOutputsUnscaled);
        // just compute
        fval = computeNegLogLikelihood(fun->getAllMeasurements(), modelOutputsUnscaled);
    }

    return functionEvaluationSuccess;
}

double HierachicalOptimizationWrapper::computeNegLogLikelihood(std::vector <std::vector<double>> const& measurements,
                                                               const std::vector<std::vector<double> > &modelOutputsScaled) const {
    double llh = 0.0;

    for (int conditionIdx = 0; conditionIdx < numConditions; ++conditionIdx) {
        llh += computeNegLogLikelihood(measurements[conditionIdx], modelOutputsScaled[conditionIdx]);
    }

    return llh;
}
double HierachicalOptimizationWrapper::computeNegLogLikelihood(std::vector<double> const& measurements,
                                                               std::vector<double> const& modelOutputsScaled) const {
    double llh = 0.0;
    constexpr double pi = atan(1)*4;

    double sigmaSquared = 1.0; // NOTE: no user-JobResultAmiciSimulationprovided sigma supported at the moment

    // No need to pay respect to timepoints and observables, as long as the order is the same for both measurements and simulations
    RELEASE_ASSERT(measurements.size() == modelOutputsScaled.size(), "measurement/simulation output dimension mismatch");
    for(int i = 0; (unsigned) i < measurements.size(); ++i) {
        double mes = measurements[i];
        if(!std::isnan(mes)) {
            double sim = modelOutputsScaled[i];
            assert(!std::isnan(sim));
            double diff = mes - sim;
            diff *= diff;
            llh += log(2.0 * pi * sigmaSquared) + diff / sigmaSquared;

        }
    }

    llh /= 2.0;
    return llh;
}


std::vector<double> HierachicalOptimizationWrapper::spliceParameters(const double * const reducedParameters, int numReduced, const std::vector<double> &scalingFactors, const std::vector<double> &offsetParameters) const {
    std::vector<double> fullParameters(numReduced + scalingFactors.size() + offsetParameters.size());
    int idxScaling = 0;
    int idxOffset = 0;
    int idxRegular = 0;

    for(int i = 0; i < (signed) fullParameters.size(); ++i) {
        if((unsigned)idxScaling < proportionalityFactorIndices.size() && proportionalityFactorIndices[idxScaling] == i)
            fullParameters[i] = scalingFactors[idxScaling++];
        else if((unsigned)idxOffset < offsetParameterIndices.size() && offsetParameterIndices[idxOffset] == i)
            fullParameters[i] = offsetParameters[idxOffset++];
        else if(idxRegular < numReduced)
            fullParameters[i] = reducedParameters[idxRegular++];
        else
            throw std::exception();
    }

    return fullParameters;
}

int HierachicalOptimizationWrapper::numParameters() const {
    return fun->numParameters() - numScalingFactors() - numOffsetParameters();
}

int HierachicalOptimizationWrapper::numScalingFactors() const {
    return proportionalityFactorIndices.size();
}

const std::vector<int> &HierachicalOptimizationWrapper::getProportionalityFactorIndices() const
{
    return proportionalityFactorIndices;
}

int HierachicalOptimizationWrapper::numOffsetParameters() const {
    return offsetParameterIndices.size();
}

const std::vector<int> &HierachicalOptimizationWrapper::getOffsetParameterIndices() const
{
    return offsetParameterIndices;
}

//
// merge with OffsetParameterHdf5Reader?
ScalingFactorHdf5Reader::ScalingFactorHdf5Reader(H5::H5File file, std::string rootPath)
    :file(file), rootPath(rootPath)
{
    scalingParameterIndicesPath = rootPath + "/scalingParameterIndices";
    mapPath = rootPath + "/scalingParametersMapToObservables";
    readFromFile();
}

std::vector<int> ScalingFactorHdf5Reader::getConditionsForScalingParameter(int scalingIdx) const {
    std::vector<int> result;
    result.reserve(mapping[scalingIdx].size());
    for (auto const& kvp : mapping[scalingIdx])
        result.push_back(kvp.first);
    return result;
}



const std::vector<int> &ScalingFactorHdf5Reader::getObservablesForScalingParameter(int scalingIdx, int conditionIdx) const {
    return mapping[scalingIdx].at(conditionIdx);
}

std::vector<int> ScalingFactorHdf5Reader::getProportionalityFactorIndices() {
    auto lock = hdf5MutexGetLock();
    std::vector<int> proportionalityFactorIndices;
    H5_SAVE_ERROR_HANDLER; // don't show error if dataset is missing
    try {
        auto dataset = file.openDataSet(scalingParameterIndicesPath);
        auto dataspace = dataset.getSpace();

        auto ndims = dataspace.getSimpleExtentNdims();
        if(ndims != 1)
            throw ParPEException("Invalid dimension in getProportionalityFactorIndices.");
        hsize_t numScalings = 0;
        dataspace.getSimpleExtentDims(&numScalings);

        proportionalityFactorIndices.resize(numScalings);
        dataset.read(proportionalityFactorIndices.data(), H5::PredType::NATIVE_INT);
    } catch (H5::FileIException e) {
    }
    H5_RESTORE_ERROR_HANDLER;

    return proportionalityFactorIndices;
}

void ScalingFactorHdf5Reader::readFromFile() {
    auto lock = hdf5MutexGetLock();
    H5_SAVE_ERROR_HANDLER;
    try {
        auto dataset = file.openDataSet(mapPath);
        auto dataspace = dataset.getSpace();

        auto attribute = dataset.openAttribute("numScalings");
        auto type = attribute.getDataType();
        int numScalings = 0;
        attribute.read(type, &numScalings);
        if(numScalings == 0)
            return;

        auto ndims = dataspace.getSimpleExtentNdims();
        if(ndims != 2)
            throw ParPEException("Invalid dimension for scaling parameter map, expected 2.");
        hsize_t dims[ndims];
        dataspace.getSimpleExtentDims(dims);
        hsize_t const nRows = dims[0];
        hsize_t const nCols = dims[1];
        if(nRows && nCols != 3)
            throw ParPEException("Invalid dimension for scaling parameter map, expected 2.");

        // column indices in dataspace
        constexpr int scalingCol = 0;
        constexpr int conditionCol = 1;
        constexpr int observableCol = 2;

        std::vector<int> rawMap(nRows * nCols);
        dataset.read(rawMap.data(), H5::PredType::NATIVE_INT);

        mapping.resize(numScalings);
        for(int i = 0; (unsigned)i < nRows; ++i) {
            int scalingIdx = rawMap[i * nCols + scalingCol];
            int conditionIdx = rawMap[i * nCols + conditionCol];
            int observableIdx = rawMap[i * nCols + observableCol];
            mapping[scalingIdx][conditionIdx].push_back(observableIdx);
        }
    } catch (H5::FileIException e) {
        return;
    }
    H5_RESTORE_ERROR_HANDLER;

}

// Offsets
OffsetParameterHdf5Reader::OffsetParameterHdf5Reader(H5::H5File file, std::string rootPath)
    :file(file), rootPath(rootPath)
{
    offsetParameterIndicesPath = rootPath + "/offsetParameterIndices";
    mapPath = rootPath + "/offsetParametersMapToObservables";
    readFromFile();
}

std::vector<int> OffsetParameterHdf5Reader::getConditionsForOffsetParameter(int offsetIdx) const {
    // TODO adapt to Offsets
    std::vector<int> result;
    result.reserve(mapping[offsetIdx].size());
    for (auto const& kvp : mapping[offsetIdx])
        result.push_back(kvp.first);
    return result;
}

const std::vector<int> &OffsetParameterHdf5Reader::getObservablesForOffsetParameter(int offsetIdx, int conditionIdx) const {
    // TODO adapt to Offsets
    return mapping[offsetIdx].at(conditionIdx);
}

std::vector<int> OffsetParameterHdf5Reader::getOffsetParameterIndices() {
    auto lock = hdf5MutexGetLock();
    std::vector<int> offsetParameterIndices;
    H5_SAVE_ERROR_HANDLER; // don't show error if dataset is missing
    try {
        auto dataset = file.openDataSet(offsetParameterIndicesPath);
        auto dataspace = dataset.getSpace();

        auto ndims = dataspace.getSimpleExtentNdims();
        if(ndims != 1)
            throw ParPEException("Invalid dimension in getOffsetParameterIndices.");
        hsize_t numOffsets = 0;
        dataspace.getSimpleExtentDims(&numOffsets);

        offsetParameterIndices.resize(numOffsets);
        dataset.read(offsetParameterIndices.data(), H5::PredType::NATIVE_INT);
    } catch (H5::FileIException e) {
    }
    H5_RESTORE_ERROR_HANDLER;

    return offsetParameterIndices;
}

void OffsetParameterHdf5Reader::readFromFile() {
    auto lock = hdf5MutexGetLock();
    H5_SAVE_ERROR_HANDLER;
    try {
        auto dataset = file.openDataSet(mapPath);
        auto dataspace = dataset.getSpace();

        auto attribute = dataset.openAttribute("numOffsets");
        auto type = attribute.getDataType();
        int numOffsets = 0;
        attribute.read(type, &numOffsets);
        if(numOffsets == 0)
            return;

        auto ndims = dataspace.getSimpleExtentNdims();
        if(ndims != 2)
            throw ParPEException("Invalid dimension for offset parameter map, expected 2.");
        hsize_t dims[ndims];
        dataspace.getSimpleExtentDims(dims);
        hsize_t const nRows = dims[0];
        hsize_t const nCols = dims[1];
        if(nRows && nCols != 3)
            throw ParPEException("Invalid dimension for offset parameter map, expected 2.");

        // column indices in dataspace
        constexpr int offsetCol = 0;
        constexpr int conditionCol = 1;
        constexpr int observableCol = 2;

        std::vector<int> rawMap(nRows * nCols);
        dataset.read(rawMap.data(), H5::PredType::NATIVE_INT);

        mapping.resize(numOffsets);
        for(int i = 0; (unsigned)i < nRows; ++i) {
            int offsetIdx = rawMap[i * nCols + offsetCol];
            int conditionIdx = rawMap[i * nCols + conditionCol];
            int observableIdx = rawMap[i * nCols + observableCol];
            mapping[offsetIdx][conditionIdx].push_back(observableIdx);
        }
    } catch (H5::FileIException e) {
        return;
    }
    H5_RESTORE_ERROR_HANDLER;

}

//
HierachicalOptimizationProblemWrapper::HierachicalOptimizationProblemWrapper(
        std::unique_ptr<OptimizationProblem> problemToWrap,
        const MultiConditionDataProvider *dataProvider)
    : wrappedProblem(std::move(problemToWrap))
{
    auto wrappedFun = dynamic_cast<SummedGradientFunctionGradientFunctionAdapter<int>*>(wrappedProblem->costFun.get());

    costFun.reset(new HierachicalOptimizationWrapper(
                      std::unique_ptr<AmiciSummedGradientFunction<int>>(
                          dynamic_cast<AmiciSummedGradientFunction<int>*>(wrappedFun->gradFun.get())),
                      std::make_unique<ScalingFactorHdf5Reader>(dataProvider->file, "/"),
                      std::make_unique<OffsetParameterHdf5Reader>(dataProvider->file, "/"),
                      dataProvider->getNumberOfConditions(),
                      dataProvider->model->nytrue,
                      dataProvider->model->nt(), ParameterTransformation::log10,
                      ErrorModel::normal));
}

HierachicalOptimizationProblemWrapper::~HierachicalOptimizationProblemWrapper()
{
    // Avoid double delete. This will be destroyed when wrappedProblem goes out of scope!
    dynamic_cast<HierachicalOptimizationWrapper *>(costFun.get())->fun.release();
}

void HierachicalOptimizationProblemWrapper::fillInitialParameters(double *buffer) const
{
    std::vector<double> full(wrappedProblem->costFun->numParameters());
    wrappedProblem->fillInitialParameters(full.data());
    fillFilteredParams(full, buffer);
}

void HierachicalOptimizationProblemWrapper::fillParametersMax(double *buffer) const
{
    std::vector<double> full(wrappedProblem->costFun->numParameters());
    wrappedProblem->fillParametersMax(full.data());
    fillFilteredParams(full, buffer);
}

void HierachicalOptimizationProblemWrapper::fillParametersMin(double *buffer) const
{
    std::vector<double> full(wrappedProblem->costFun->numParameters());
    wrappedProblem->fillParametersMin(full.data());
    fillFilteredParams(full, buffer);
}

void HierachicalOptimizationProblemWrapper::fillFilteredParams(const std::vector<double> &fullParams, double *buffer) const
{
    auto hierarchical = dynamic_cast<HierachicalOptimizationWrapper *>(costFun.get());
    auto proportionalityFactorIndices = hierarchical->getProportionalityFactorIndices();
    auto offsetParameterIndices = hierarchical->getOffsetParameterIndices();
    auto scalingParameterIndices = proportionalityFactorIndices;
    scalingParameterIndices.insert(scalingParameterIndices.end(),offsetParameterIndices.begin(),offsetParameterIndices.end());
    std::sort(scalingParameterIndices.begin(),scalingParameterIndices.end());
    parpe::fillFilteredParams(fullParams, scalingParameterIndices, buffer);
}

void fillFilteredParams(const std::vector<double> &fullParams, std::vector<int> const& sortedFilterIndices, double *buffer)
{
    // adapt to offsets
    unsigned int nextFilterIdx = 0;
    unsigned int bufferIdx = 0;
    for(int i = 0; (unsigned)i < fullParams.size(); ++i) {
        if(nextFilterIdx < sortedFilterIndices.size()
                && sortedFilterIndices[nextFilterIdx] == i) {
            // skip
            ++nextFilterIdx;
        } else {
            // copy
            buffer[bufferIdx] = fullParams[i];
            ++bufferIdx;
        }
    }
    assert(nextFilterIdx == sortedFilterIndices.size());
    assert(bufferIdx == (unsigned) fullParams.size() - sortedFilterIndices.size());

}


} // namespace parpe
