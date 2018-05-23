#include "hierarchicalOptimization.h"

#include <assert.h>
#include <exception>

#ifdef __INTEL_COMPILER
// constexpr did not work on icc (ICC) 16.0.4 20160811
#define constexpr
#endif

namespace parpe {


HierachicalOptimizationWrapper::HierachicalOptimizationWrapper(std::unique_ptr<AmiciSummedGradientFunction<int> > fun,
                                                               int numConditions, int numObservables, int numTimepoints)
    : fun(std::move(fun)),
      numConditions(numConditions),
      numObservables(numObservables),
      numTimepoints(numTimepoints)
{
    scalingReader = std::make_unique<AnalyticalParameterHdf5Reader>();
    offsetReader = std::make_unique<AnalyticalParameterHdf5Reader>();
    init();
}

HierachicalOptimizationWrapper::HierachicalOptimizationWrapper(std::unique_ptr<AmiciSummedGradientFunction<int> > fun,
                                                               H5::H5File const& file, std::string hdf5RootPath,
                                                               int numConditions, int numObservables, int numTimepoints,
                                                               ErrorModel errorModel)
    : fun(std::move(fun)),
      numConditions(numConditions),
      numObservables(numObservables),
      numTimepoints(numTimepoints),
      errorModel(errorModel)
{
    scalingReader = std::make_unique<AnalyticalParameterHdf5Reader>(file,
                                                                    hdf5RootPath + "/scalingParameterIndices",
                                                                    hdf5RootPath + "/scalingParametersMapToObservables");
    offsetReader = std::make_unique<AnalyticalParameterHdf5Reader>(file,
                                                                   hdf5RootPath + "/offsetParameterIndices",
                                                                   hdf5RootPath + "/offsetParametersMapToObservables");

    init();
}

HierachicalOptimizationWrapper::HierachicalOptimizationWrapper(
        std::unique_ptr<AmiciSummedGradientFunction<int> > fun,
        std::unique_ptr<parpe::AnalyticalParameterProvider> scalingReader,
        std::unique_ptr<parpe::AnalyticalParameterProvider> offsetReader,
        int numConditions, int numObservables, int numTimepoints,
        ErrorModel errorModel)
    : fun(std::move(fun)),
      scalingReader(std::move(scalingReader)),
      offsetReader(std::move(offsetReader)),
      numConditions(numConditions),
      numObservables(numObservables),
      numTimepoints(numTimepoints),
      errorModel(errorModel)
{
    init();
}


void HierachicalOptimizationWrapper::init() {
    if(errorModel != ErrorModel::normal) {
        throw ParPEException("Only gaussian noise is supported so far.");
    }

    // some functions currently expect these lists to be sorted, therefore ensure sorting right away
    // (if sorting here, also need to reorder/reindex scalingFactorIdx in mapping table -> difficult)
    proportionalityFactorIndices = this->scalingReader->getOptimizationParameterIndices();
    RELEASE_ASSERT(std::is_sorted(this->proportionalityFactorIndices.begin(),
                                  this->proportionalityFactorIndices.end()), "");

    offsetParameterIndices = this->offsetReader->getOptimizationParameterIndices();
    RELEASE_ASSERT(std::is_sorted(this->offsetParameterIndices.begin(),
                                  this->offsetParameterIndices.end()), "");
}


FunctionEvaluationStatus HierachicalOptimizationWrapper::evaluate(const double * const parameters,
                                                                  double &fval,
                                                                  double *gradient) const {
    if(numScalingFactors() == 0 && numOffsetParameters() == 0) {
        // nothing to do, just pass through

        // evaluate for all conditions
        std::vector<int> dataIndices(numConditions);
        std::iota(dataIndices.begin(), dataIndices.end(), 0);

        return fun->evaluate(parameters, dataIndices, fval, gradient);
    }

    // evaluate with scaling parameters set to 1 and offsets to 0
    auto modelOutput = getUnscaledModelOutputs(parameters);

    auto measurements = fun->getAllMeasurements();

    // compute correct scaling factors analytically
    auto scalings = computeAnalyticalScalings(measurements, modelOutput);

    // compute correct offset parameters analytically
    auto offsets = computeAnalyticalOffsets(measurements, modelOutput);
    // std::cout << "offsets:" << offsets;

    // evaluate with analytical scaling parameters
    auto status = evaluateWithOptimalParameters(parameters, scalings, offsets,
                                                measurements, modelOutput,
                                                fval, gradient);

    return status;
}

std::vector<double> HierachicalOptimizationWrapper::getDefaultScalingFactors() const
{
    auto result = std::vector<double>(numScalingFactors());

    for(int i = 0; i < numScalingFactors(); ++i) {
        result[i] = getDefaultScalingFactor(
                    fun->getParameterScaling(proportionalityFactorIndices[i]));
    }

    return result;
}

std::vector<double> HierachicalOptimizationWrapper::getDefaultOffsetParameters() const
{
    auto result = std::vector<double>(numOffsetParameters());

    for(int i = 0; i < numOffsetParameters(); ++i) {
        result[i] = getDefaultOffsetParameter(
                    fun->getParameterScaling(offsetParameterIndices[i]));
    }

    return result;
}

std::vector<std::vector<double> > HierachicalOptimizationWrapper::getUnscaledModelOutputs(const double * const reducedParameters) const {
    // run simulations, collect outputs
    auto scalingDummy = getDefaultScalingFactors();
    auto offsetDummy = getDefaultOffsetParameters();

    // splice hidden scaling parameter and external parameters
    auto fullParameters = spliceParameters(reducedParameters, numParameters(),
                                           proportionalityFactorIndices, offsetParameterIndices,
                                           scalingDummy, offsetDummy);

    std::vector<std::vector<double> > modelOutput(numConditions);
    fun->getModelOutputs(fullParameters.data(), modelOutput);

    return modelOutput;
}

std::vector<double> HierachicalOptimizationWrapper::computeAnalyticalScalings(
        std::vector<std::vector<double>> const& measurements,
        std::vector<std::vector<double>> &modelOutputsUnscaled) const
{
    // NOTE: does not handle replicates, assumes normal distribution, does not compute sigmas

    int numProportionalityFactors = proportionalityFactorIndices.size();
    std::vector<double> proportionalityFactors(numProportionalityFactors);

    for(int scalingIdx = 0; scalingIdx < numProportionalityFactors; ++scalingIdx) {
        proportionalityFactors[scalingIdx] = parpe::computeAnalyticalScalings(scalingIdx,
                                                                              fun->getParameterScaling(proportionalityFactorIndices[scalingIdx]),
                                                                              modelOutputsUnscaled, measurements,
                                                                              *scalingReader, numObservables, numTimepoints);
    }

    return proportionalityFactors;
}

void HierachicalOptimizationWrapper::applyOptimalScalings(std::vector<double> const& proportionalityFactors,
                                                          std::vector<std::vector<double> > &modelOutputs) const {

    for(int i = 0; (unsigned) i < proportionalityFactors.size(); ++i) {
        double scaling = getUnscaledParameter(proportionalityFactors[i],
                                              fun->getParameterScaling(proportionalityFactorIndices[i]));

        applyOptimalScaling(i, scaling, modelOutputs,
                            *scalingReader, numObservables, numTimepoints);
    }
}


std::vector<double> HierachicalOptimizationWrapper::computeAnalyticalOffsets(std::vector<std::vector<double>> const& measurements,
                                                                             std::vector<std::vector<double> > &modelOutputsUnscaled) const {
    // NOTE: does not handle replicates, assumes normal distribution, does not compute sigmas

    int numOffsetParameters = offsetParameterIndices.size();
    std::vector<double> offsetParameters(numOffsetParameters);

    for(int i = 0; i < numOffsetParameters; ++i) {
        parpe::computeAnalyticalOffsets(i,
                                        fun->getParameterScaling(offsetParameterIndices[i]),
                                        modelOutputsUnscaled, measurements,
                                        *offsetReader, numObservables, numTimepoints);
    }

    return offsetParameters;
}

void HierachicalOptimizationWrapper::applyOptimalOffsets(std::vector<double> const& offsetParameters,
                                                         std::vector<std::vector<double> > &modelOutputs) const {

    // TODO: do we need to distinguish here?
    for(int i = 0; (unsigned) i < offsetParameters.size(); ++i) {
        double offset = getUnscaledParameter(offsetParameters[i],
                                             fun->getParameterScaling(offsetParameterIndices[i]));
        applyOptimalOffset(i, offset, modelOutputs,
                           *offsetReader, numObservables, numTimepoints);
    }
}



FunctionEvaluationStatus HierachicalOptimizationWrapper::evaluateWithOptimalParameters(
        // adapt to offsets
        const double * const reducedParameters,
        const std::vector<double> &scalings,
        const std::vector<double> &offsets,
        std::vector<std::vector<double>> const& measurements,
        std::vector<std::vector<double> > &modelOutputsUnscaled,
        double &fval, double *gradient) const {

    if(gradient) {
        // simulate with updated theta for sensitivities
        // splice parameter vector we get from optimizer with analytically computed parameters
        auto fullParameters = spliceParameters(reducedParameters, numParameters(),
                                               proportionalityFactorIndices, offsetParameterIndices,
                                               scalings, offsets);

        // simualte all datasets
        std::vector<int> dataIndices(numConditions);
        std::iota(dataIndices.begin(), dataIndices.end(), 0);

        // Need intermediary buffer because optimizer expects fewer parameters than `fun` delivers
        std::vector<double> fullGradient(fullParameters.size());
        auto status = fun->evaluate(fullParameters.data(), dataIndices, fval, fullGradient.data());
        if(status != functionEvaluationSuccess)
            return status;

        // Filter gradient for those parameters expected by the optimizer
        auto analyticalParameterIndices = getAnalyticalParameterIndices();
        fillFilteredParams(fullGradient, analyticalParameterIndices, gradient);
    } else {
        // Scale model outputs
        applyOptimalScalings(scalings, modelOutputsUnscaled);
        applyOptimalOffsets(offsets, modelOutputsUnscaled);
        // ... to compute negative log-likelihood
        fval = computeNegLogLikelihood(measurements, modelOutputsUnscaled);
    }

    return functionEvaluationSuccess;
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


AnalyticalParameterHdf5Reader::AnalyticalParameterHdf5Reader(H5::H5File const& file,
                                                             std::string scalingParameterIndicesPath,
                                                             std::string mapPath)
    : mapPath(mapPath),
      analyticalParameterIndicesPath(scalingParameterIndicesPath)
{
    auto lock = hdf5MutexGetLock();
    this->file = file; // copy while mutex is locked!
    readParameterConditionObservableMappingFromFile();
}


int HierachicalOptimizationWrapper::numOffsetParameters() const {
    return offsetParameterIndices.size();
}

const std::vector<int> &HierachicalOptimizationWrapper::getOffsetParameterIndices() const
{
    return offsetParameterIndices;
}

std::vector<int> HierachicalOptimizationWrapper::getAnalyticalParameterIndices() const
{
    auto combinedIndices = proportionalityFactorIndices;
    combinedIndices.insert(combinedIndices.end(), offsetParameterIndices.begin(), offsetParameterIndices.end());
    std::sort(combinedIndices.begin(), combinedIndices.end());

    return combinedIndices;
}

std::vector<int> AnalyticalParameterHdf5Reader::getConditionsForParameter(int parameterIndex) const {
    std::vector<int> result;
    result.reserve(mapping[parameterIndex].size());
    for (auto const& kvp : mapping[parameterIndex])
        result.push_back(kvp.first);
    return result;
}

const std::vector<int> &AnalyticalParameterHdf5Reader::getObservablesForParameter(
        int parameterIndex, int conditionIdx) const {
    return mapping[parameterIndex].at(conditionIdx);
}


std::vector<int> AnalyticalParameterHdf5Reader::getOptimizationParameterIndices() const {
    auto lock = hdf5MutexGetLock();
    std::vector<int> analyticalParameterIndices;
    H5_SAVE_ERROR_HANDLER; // don't show error if dataset is missing
    try {
        auto dataset = file.openDataSet(analyticalParameterIndicesPath);
        auto dataspace = dataset.getSpace();

        auto ndims = dataspace.getSimpleExtentNdims();
        if(ndims != 1)
            throw ParPEException("Invalid dimension in getOptimizationParameterIndices.");
        hsize_t numScalings = 0;
        dataspace.getSimpleExtentDims(&numScalings);

        analyticalParameterIndices.resize(numScalings);
        dataset.read(analyticalParameterIndices.data(), H5::PredType::NATIVE_INT);
    } catch (H5::FileIException e) {
        // we just return an empty list
    }
    H5_RESTORE_ERROR_HANDLER;

    return analyticalParameterIndices;
}

int AnalyticalParameterHdf5Reader::getNumAnalyticalParameters(H5::DataSet& dataset) const
{
    hsize_t numAnalyticalParameters = 0;
    auto lock = hdf5MutexGetLock();

    H5_SAVE_ERROR_HANDLER; // don't show error if dataset is missing
    try {
        auto dataset = file.openDataSet(analyticalParameterIndicesPath);
        auto dataspace = dataset.getSpace();
        auto ndims = dataspace.getSimpleExtentNdims();
        if(ndims != 1)
            throw ParPEException("Invalid dimension in getOptimizationParameterIndices.");
        dataspace.getSimpleExtentDims(&numAnalyticalParameters);
    } catch (H5::FileIException e) {
        // 0
    }
    H5_RESTORE_ERROR_HANDLER;

    return numAnalyticalParameters;
}

void AnalyticalParameterHdf5Reader::readParameterConditionObservableMappingFromFile() {
    auto lock = hdf5MutexGetLock();
    H5_SAVE_ERROR_HANDLER;
    try {
        auto dataset = file.openDataSet(mapPath);
        int numScalings = getNumAnalyticalParameters(dataset);
        if(numScalings == 0)
            return;

        // column indices in dataspace
        constexpr int parameterCol = 0;
        constexpr int conditionCol = 1;
        constexpr int observableCol = 2;

        mapping.resize(numScalings);

        hsize_t nRows = 0, nCols = 0;
        auto rawMap = readRawMap(dataset, nRows, nCols);

        for(int i = 0; (unsigned)i < nRows; ++i) {
            int scalingIdx = rawMap[i * nCols + parameterCol];
            int conditionIdx = rawMap[i * nCols + conditionCol];
            int observableIdx = rawMap[i * nCols + observableCol];
            mapping[scalingIdx][conditionIdx].push_back(observableIdx);
        }
    } catch (H5::FileIException e) {
        return;
    }
    H5_RESTORE_ERROR_HANDLER;

}

std::vector<int> AnalyticalParameterHdf5Reader::readRawMap(H5::DataSet& dataset, hsize_t& nRows, hsize_t& nCols)
{
    auto dataspace = dataset.getSpace();
    auto ndims = dataspace.getSimpleExtentNdims();
    if(ndims != 2)
        throw ParPEException("Invalid dimension for analytical parameter map, expected 2.");
    hsize_t dims[ndims];
    dataspace.getSimpleExtentDims(dims);
    nRows = dims[0];
    nCols = dims[1];
    if(nRows && nCols != 3)
        throw ParPEException("Invalid dimension for analytical parameter map, expected 2.");

    std::vector<int> rawMap(nRows * nCols);
    dataset.read(rawMap.data(), H5::PredType::NATIVE_INT);

    return rawMap;
}

//
HierachicalOptimizationProblemWrapper::HierachicalOptimizationProblemWrapper(
        std::unique_ptr<OptimizationProblem> problemToWrap,
        const MultiConditionDataProviderHDF5 *dataProvider)
    : wrappedProblem(std::move(problemToWrap))
{
    auto wrappedFun = dynamic_cast<SummedGradientFunctionGradientFunctionAdapter<int>*>(
                wrappedProblem->costFun.get());

    auto model = dataProvider->getModel();

    costFun.reset(new HierachicalOptimizationWrapper(
                      std::unique_ptr<AmiciSummedGradientFunction<int>>(
                          dynamic_cast<AmiciSummedGradientFunction<int>*>(wrappedFun->getWrappedFunction())),
                      dataProvider->getHdf5FileId(), "/",
                      dataProvider->getNumberOfConditions(),
                      model->nytrue,
                      model->nt(),
                      ErrorModel::normal));
}

HierachicalOptimizationProblemWrapper::HierachicalOptimizationProblemWrapper(std::unique_ptr<OptimizationProblem> problemToWrap,
                                                                             std::unique_ptr<HierachicalOptimizationWrapper> costFun)
    : OptimizationProblem(std::move(costFun)),
      wrappedProblem(std::move(problemToWrap))
{

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
    auto combinedIndices = hierarchical->getAnalyticalParameterIndices();
    parpe::fillFilteredParams(fullParams, combinedIndices, buffer);
}

void fillFilteredParams(std::vector<double> const& valuesToFilter,
                        std::vector<int> const& sortedIndicesToExclude,
                        double *result)
{
    // adapt to offsets
    unsigned int nextFilterIdx = 0;
    unsigned int resultIdx = 0;
    for(int i = 0; (unsigned)i < valuesToFilter.size(); ++i) {
        if(nextFilterIdx < sortedIndicesToExclude.size()
                && sortedIndicesToExclude[nextFilterIdx] == i) {
            // skip
            ++nextFilterIdx;
        } else {
            // copy
            result[resultIdx] = valuesToFilter[i];
            ++resultIdx;
        }
    }
    RELEASE_ASSERT(nextFilterIdx == sortedIndicesToExclude.size(), "");
    RELEASE_ASSERT(resultIdx == (unsigned) valuesToFilter.size() - sortedIndicesToExclude.size(), "");
}

double getDefaultScalingFactor(amici::AMICI_parameter_scaling scaling)
{
    switch (scaling) {
    case amici::AMICI_SCALING_NONE:
        return 1.0;
    case amici::AMICI_SCALING_LOG10:
        return 0.0;
    default:
        throw ParPEException("Parameter scaling must be AMICI_SCALING_LOG10 or AMICI_SCALING_NONE.");
    }
}

double getDefaultOffsetParameter(amici::AMICI_parameter_scaling scaling)
{
    switch (scaling) {
    case amici::AMICI_SCALING_NONE:
        return 0.0;
    case amici::AMICI_SCALING_LOG10:
        return -INFINITY;
    default:
        throw ParPEException("Parameter scaling must be AMICI_SCALING_LOG10 or AMICI_SCALING_NONE.");
    }
}


double computeAnalyticalScalings(int scalingIdx, amici::AMICI_parameter_scaling scale,
                                 const std::vector<std::vector<double> > &modelOutputsUnscaled,
                                 const std::vector<std::vector<double> > &measurements,
                                 AnalyticalParameterProvider& scalingReader,
                                 int numObservables, int numTimepoints) {

    auto dependentConditions = scalingReader.getConditionsForParameter(scalingIdx);

    double enumerator = 0.0;
    double denominator = 0.0;

    for (auto const conditionIdx: dependentConditions) {
        auto dependentObservables = scalingReader.getObservablesForParameter(scalingIdx, conditionIdx);
        for(auto const observableIdx: dependentObservables) {
            for(int timeIdx = 0; timeIdx < numTimepoints; ++timeIdx) {
                double mes = measurements[conditionIdx][observableIdx + timeIdx * numObservables];
                if(!std::isnan(mes)) {
                    // NOTE: this must be in sync with data ordering in AMICI (assumes row-major)
                    double sim = modelOutputsUnscaled[conditionIdx][observableIdx + timeIdx * numObservables];
                    assert(!std::isnan(sim));
                    enumerator += sim * mes;
                    denominator += sim * sim;
                }
            }
        }
    }

    double proportionalityFactor = getScaledParameter(enumerator / denominator, scale);

    return proportionalityFactor;
}


double computeAnalyticalOffsets(int offsetIdx,
                                amici::AMICI_parameter_scaling scale,
                                const std::vector<std::vector<double> > &modelOutputsUnscaled,
                                const std::vector<std::vector<double> > &measurements,
                                AnalyticalParameterProvider& offsetReader,
                                int numObservables, int numTimepoints) {
    auto dependentConditions = offsetReader.getConditionsForParameter(offsetIdx);

    double enumerator = 0.0;
    double denominator = 0.0;

    for (auto const conditionIdx: dependentConditions) {
        auto dependentObservables = offsetReader.getObservablesForParameter(offsetIdx, conditionIdx);
        for(auto const observableIdx: dependentObservables) {
            for(int timeIdx = 0; timeIdx < numTimepoints; ++timeIdx) {
                double mes = measurements[conditionIdx][observableIdx + timeIdx * numObservables];
                if(!std::isnan(mes)) {
                    double sim = modelOutputsUnscaled[conditionIdx][observableIdx + timeIdx * numObservables];
                    assert(!std::isnan(sim));
                    enumerator += mes - sim;
                    denominator += 1.0;
                }
            }
        }
    }

    double offsetParameter = getScaledParameter(enumerator / denominator, scale);

    // TODO ensure positivity!
    return offsetParameter;
}


void applyOptimalScaling(int scalingIdx, double scalingLin,
                         std::vector<std::vector<double> > &modelOutputs,
                         AnalyticalParameterProvider& scalingReader,
                         int numObservables, int numTimepoints) {
    auto dependentConditions = scalingReader.getConditionsForParameter(scalingIdx);
    for (auto const conditionIdx: dependentConditions) {
        auto dependentObservables = scalingReader.getObservablesForParameter(scalingIdx, conditionIdx);
        for(auto const observableIdx: dependentObservables) {
            assert(observableIdx < numObservables);
            for(int timeIdx = 0; timeIdx < numTimepoints; ++timeIdx) {
                // NOTE: this must be in sync with data ordering in AMICI (assumes row-major)
                modelOutputs[conditionIdx][observableIdx + timeIdx * numObservables] *= scalingLin;
            }
        }
    }
}

double getScaledParameter(double parameter, amici::AMICI_parameter_scaling scale) {
    switch (scale) {
    case amici::AMICI_SCALING_NONE:
        return parameter;
    case amici::AMICI_SCALING_LOG10:
        return log10(parameter);
    default:
        throw(ParPEException("Parameter scaling must be AMICI_SCALING_LOG10 or AMICI_SCALING_NONE."));
    }
}

double getUnscaledParameter(double parameter, amici::AMICI_parameter_scaling scale) {
    switch (scale) {
    case amici::AMICI_SCALING_NONE:
        return parameter;
    case amici::AMICI_SCALING_LOG10:
        return pow(10, parameter);
    default:
        throw ParPEException("Parameter scaling must be AMICI_SCALING_LOG10 or AMICI_SCALING_NONE.");
    }
}


void applyOptimalOffset(int offsetIdx, double offsetLin,
                        std::vector<std::vector<double> > &modelOutputs,
                        AnalyticalParameterProvider& offsetReader,
                        int numObservables, int numTimepoints) {
    auto dependentConditions = offsetReader.getConditionsForParameter(offsetIdx);
    for (auto const conditionIdx: dependentConditions) {
        auto dependentObservables = offsetReader.getObservablesForParameter(offsetIdx, conditionIdx);
        for(auto const observableIdx: dependentObservables) {
            assert(observableIdx < numObservables);
            for(int timeIdx = 0; timeIdx < numTimepoints; ++timeIdx) {
                modelOutputs[conditionIdx][observableIdx + timeIdx * numObservables] += offsetLin;
            }
        }
    }
}



std::vector<double> spliceParameters(const double * const reducedParameters, int numReduced,
                                     const std::vector<int> &proportionalityFactorIndices,
                                     const std::vector<int> &offsetParameterIndices,
                                     const std::vector<double> &scalingFactors,
                                     const std::vector<double> &offsetParameters) {
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


double computeNegLogLikelihood(std::vector <std::vector<double>> const& measurements,
                               const std::vector<std::vector<double> > &modelOutputsScaled) {
    RELEASE_ASSERT(measurements.size() == modelOutputsScaled.size(), "");

    double nllh = 0.0;

    for (int conditionIdx = 0; (unsigned) conditionIdx < measurements.size(); ++conditionIdx) {
        nllh += computeNegLogLikelihood(measurements[conditionIdx], modelOutputsScaled[conditionIdx]);
    }

    return nllh;
}

double computeNegLogLikelihood(std::vector<double> const& measurements,
                               std::vector<double> const& modelOutputsScaled) {
    double nllh = 0.0;
    constexpr double pi = atan(1)*4.0;

    double sigmaSquared = 1.0; // NOTE: TODO: no user-JobResultAmiciSimulationprovided sigma supported at the moment

    // No need to pay respect to timepoints and observables, as long as the order is the same for both measurements and simulations
    RELEASE_ASSERT(measurements.size() == modelOutputsScaled.size(), "measurement/simulation output dimension mismatch");
    for(int i = 0; (unsigned) i < measurements.size(); ++i) {
        double mes = measurements[i];
        if(!std::isnan(mes)) {
            double sim = modelOutputsScaled[i];
            assert(!std::isnan(sim));
            double diff = mes - sim;
            diff *= diff;
            nllh += log(2.0 * pi * sigmaSquared) + diff / sigmaSquared;
        }
    }

    nllh /= 2.0;
    return nllh;
}

std::vector<int> AnalyticalParameterProviderDefault::getConditionsForParameter(int parameterIndex) const {
    return conditionsForParameter[parameterIndex];
}

const std::vector<int> &AnalyticalParameterProviderDefault::getObservablesForParameter(int parameterIndex, int conditionIdx) const {
    return mapping[parameterIndex].at(conditionIdx);
}

std::vector<int> AnalyticalParameterProviderDefault::getOptimizationParameterIndices() const {
    return optimizationParameterIndices;
}


} // namespace parpe
