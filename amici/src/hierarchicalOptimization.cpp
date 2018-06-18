#include "hierarchicalOptimization.h"
#include <logging.h>

#include <assert.h>
#include <exception>

#ifdef __INTEL_COMPILER
// constexpr did not work on icc (ICC) 16.0.4 20160811
#define constexpr
#endif

namespace parpe {


HierachicalOptimizationWrapper::HierachicalOptimizationWrapper(
        std::unique_ptr<AmiciSummedGradientFunction<int> > fun,
        int numConditions, int numObservables, int numTimepoints)
    : fun(std::move(fun)),
      numConditions(numConditions),
      numObservables(numObservables),
      numTimepoints(numTimepoints)
{
    scalingReader = std::make_unique<AnalyticalParameterHdf5Reader>();
    offsetReader = std::make_unique<AnalyticalParameterHdf5Reader>();
    sigmaReader = std::make_unique<AnalyticalParameterHdf5Reader>();
    init();
}

HierachicalOptimizationWrapper::HierachicalOptimizationWrapper(
        std::unique_ptr<AmiciSummedGradientFunction<int> > fun,
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
    sigmaReader = std::make_unique<AnalyticalParameterHdf5Reader>(file,
                                                                  hdf5RootPath + "/sigmaParameterIndices",
                                                                  hdf5RootPath + "/sigmaParametersMapToObservables");

    init();
}

HierachicalOptimizationWrapper::HierachicalOptimizationWrapper(
        std::unique_ptr<AmiciSummedGradientFunction<int> > fun,
        std::unique_ptr<parpe::AnalyticalParameterProvider> scalingReader,
        std::unique_ptr<parpe::AnalyticalParameterProvider> offsetReader,
        std::unique_ptr<parpe::AnalyticalParameterProvider> sigmaReader,
        int numConditions, int numObservables, int numTimepoints,
        ErrorModel errorModel)
    : fun(std::move(fun)),
      scalingReader(std::move(scalingReader)),
      offsetReader(std::move(offsetReader)),
      sigmaReader(std::move(sigmaReader)),
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

    sigmaParameterIndices = this->sigmaReader->getOptimizationParameterIndices();
    RELEASE_ASSERT(std::is_sorted(this->sigmaParameterIndices.begin(),
                                  this->sigmaParameterIndices.end()), "");

    std::cout<<"HierachicalOptimizationWrapper parameters: "
            <<fun->numParameters()<<" total, "
           <<numParameters()<< " numerical, "
          <<proportionalityFactorIndices.size()<<" proportionality, "
         <<offsetParameterIndices.size()<<" offset, "
        <<sigmaParameterIndices.size()<<" sigma\n";
}


FunctionEvaluationStatus HierachicalOptimizationWrapper::evaluate(
        gsl::span<const double> parameters,
        double &fval,
        gsl::span<double> gradient) const {

    std::vector<double> fullParameters;
    std::vector<double> fullGradient;
    return evaluate(parameters, fval, gradient, fullParameters, fullGradient);
}

FunctionEvaluationStatus HierachicalOptimizationWrapper::evaluate(
        gsl::span<const double> reducedParameters,
        double &fval,
        gsl::span<double> gradient,
        std::vector<double> &fullParameters, std::vector<double> &fullGradient) const
{
    RELEASE_ASSERT(reducedParameters.size() == (unsigned)numParameters(), "");
    RELEASE_ASSERT(gradient.size() == 0 || gradient.size() == reducedParameters.size(), "");
    if(numProportionalityFactors() == 0 && numOffsetParameters() == 0) {
        // nothing to do, just pass through

        // evaluate for all conditions
        std::vector<int> dataIndices(numConditions);
        std::iota(dataIndices.begin(), dataIndices.end(), 0);

        return fun->evaluate(reducedParameters, dataIndices, fval, gradient);
    }

    // evaluate with scaling parameters set to 1 and offsets to 0
    auto modelOutput = getUnscaledModelOutputs(reducedParameters);

    auto measurements = fun->getAllMeasurements();

    // compute correct scaling factors analytically
    auto scalings = computeAnalyticalScalings(measurements, modelOutput);

    // compute correct offset parameters analytically
    auto offsets = computeAnalyticalOffsets(measurements, modelOutput);
    // std::cout << "offsets:" << offsets;

    // Scale model outputs
    applyOptimalScalings(scalings, modelOutput);
    applyOptimalOffsets(offsets, modelOutput);

    // needs scaled outputs
    auto sigmas = computeAnalyticalSigmas(measurements, modelOutput);

    // splice parameter vector we get from optimizer with analytically computed parameters
    fullParameters = spliceParameters(reducedParameters,
                                      proportionalityFactorIndices, offsetParameterIndices, sigmaParameterIndices,
                                      scalings, offsets, sigmas);


    // evaluate with analytical scaling parameters
    auto status = evaluateWithOptimalParameters(fullParameters, sigmas,
                                                measurements, modelOutput,
                                                fval, gradient,
                                                fullGradient);

    return status;
}

std::vector<double> HierachicalOptimizationWrapper::getDefaultScalingFactors() const
{
    auto result = std::vector<double>(numProportionalityFactors());

    for(int i = 0; i < numProportionalityFactors(); ++i) {
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

std::vector<double> HierachicalOptimizationWrapper::getDefaultSigmaParameters() const
{
    auto result = std::vector<double>(numSigmaParameters());

    for(int i = 0; i < numSigmaParameters(); ++i) {
        // default sigma is the same as default scaling
        result[i] = getDefaultScalingFactor(
                    fun->getParameterScaling(sigmaParameterIndices[i]));
    }

    return result;

}

std::vector<std::vector<double> > HierachicalOptimizationWrapper::getUnscaledModelOutputs(
        const gsl::span<const double> reducedParameters) const {
    // run simulations, collect outputs
    auto scalingDummy = getDefaultScalingFactors();
    auto offsetDummy = getDefaultOffsetParameters();
    auto sigmaDummy = getDefaultSigmaParameters();

    // splice hidden scaling parameter and external parameters
    auto fullParameters = spliceParameters(reducedParameters,
                                           proportionalityFactorIndices, offsetParameterIndices, sigmaParameterIndices,
                                           scalingDummy, offsetDummy, sigmaDummy);

    std::vector<std::vector<double> > modelOutput(numConditions);
    fun->getModelOutputs(fullParameters, modelOutput);

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
    // NOTE: does not handle replicates, assumes normal distribution

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

std::vector<double> HierachicalOptimizationWrapper::computeAnalyticalSigmas(
        const std::vector<std::vector<double> > &measurements,
        std::vector<std::vector<double> > &modelOutputsScaled) const
{
    // NOTE: does not handle replicates, assumes normal distribution

    int numSigmas = sigmaParameterIndices.size();
    std::vector<double> sigmas(numSigmas);

    for(int i = 0; i < numSigmas; ++i) {
        parpe::computeAnalyticalSigmas(
                    i,
                    fun->getParameterScaling(sigmaParameterIndices[i]),
                    modelOutputsScaled, measurements,
                    *sigmaReader, numObservables, numTimepoints);
    }
    return sigmas;
}

void HierachicalOptimizationWrapper::applyOptimalOffsets(std::vector<double> const& offsetParameters,
                                                         std::vector<std::vector<double> > &modelOutputs) const {

    for(int i = 0; (unsigned) i < offsetParameters.size(); ++i) {
        double offset = getUnscaledParameter(offsetParameters[i],
                                             fun->getParameterScaling(offsetParameterIndices[i]));
        applyOptimalOffset(i, offset, modelOutputs,
                           *offsetReader, numObservables, numTimepoints);
    }
}

void HierachicalOptimizationWrapper::fillInAnalyticalSigmas(
        std::vector<std::vector<double> > &allSigmas,
        std::vector<double> const& analyticalSigmas) const
{
    for(int sigmaParameterIdx = 0; (unsigned) sigmaParameterIdx < analyticalSigmas.size(); ++sigmaParameterIdx) {
        // sigma value will be used for likelihood computation and not passed to AMICI -> unscale
        auto sigmaParameterValue = getUnscaledParameter(
                    analyticalSigmas[sigmaParameterIdx],
                    fun->getParameterScaling(sigmaParameterIndices[sigmaParameterIdx]));

        auto dependentConditions = sigmaReader->getConditionsForParameter(sigmaParameterIdx);
        for (auto const conditionIdx: dependentConditions) {
            auto dependentObservables = sigmaReader->getObservablesForParameter(sigmaParameterIdx, conditionIdx);
            for(auto const observableIdx: dependentObservables) {
                RELEASE_ASSERT(observableIdx < numObservables, "");
                for(int timeIdx = 0; timeIdx < numTimepoints; ++timeIdx) {
                    // NOTE: this must be in sync with data ordering in AMICI (assumes row-major)
                    // if this sigma was to be estimated, data must be nan
                    RELEASE_ASSERT(std::isnan(allSigmas[conditionIdx][observableIdx + timeIdx * numObservables]), "");
                    allSigmas[conditionIdx][observableIdx + timeIdx * numObservables] = sigmaParameterValue;
                }
            }
        }
    }
}



FunctionEvaluationStatus HierachicalOptimizationWrapper::evaluateWithOptimalParameters(
        std::vector<double> const& fullParameters,
        const std::vector<double> &sigmas,
        std::vector<std::vector<double>> const& measurements,
        std::vector<std::vector<double> > &modelOutputsScaled,
        double &fval, const gsl::span<double> gradient,
        std::vector<double>& fullGradient) const {

    if(gradient.size()) {
        // simulate with updated theta for sensitivities
        // simulate all datasets
        std::vector<int> dataIndices(numConditions);
        std::iota(dataIndices.begin(), dataIndices.end(), 0);

        // Need intermediary buffer because optimizer expects fewer parameters than `fun` delivers
        fullGradient.resize(fullParameters.size());
        auto status = fun->evaluate(fullParameters, dataIndices, fval, fullGradient);
        if(status != functionEvaluationSuccess)
            return status;

        // Filter gradient for those parameters expected by the optimizer
        auto analyticalParameterIndices = getAnalyticalParameterIndices();
        fillFilteredParams(fullGradient, analyticalParameterIndices, gradient);
    } else {
        auto fullSigmaMatrices = fun->getAllSigmas();
        if(sigmaParameterIndices.size()) {
            fillInAnalyticalSigmas(fullSigmaMatrices, sigmas);
        }

        // ... to compute negative log-likelihood
        fval = computeNegLogLikelihood(measurements, modelOutputsScaled, fullSigmaMatrices);
    }

    return functionEvaluationSuccess;
}


int HierachicalOptimizationWrapper::numParameters() const {
    return fun->numParameters() - numProportionalityFactors() - numOffsetParameters() - numSigmaParameters();
}

int HierachicalOptimizationWrapper::numProportionalityFactors() const {
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

int HierachicalOptimizationWrapper::numSigmaParameters() const
{
    return sigmaParameterIndices.size();
}

const std::vector<int> &HierachicalOptimizationWrapper::getOffsetParameterIndices() const
{
    return offsetParameterIndices;
}

const std::vector<int> &HierachicalOptimizationWrapper::getSigmaParameterIndices() const
{
    return sigmaParameterIndices;
}


std::vector<int> HierachicalOptimizationWrapper::getAnalyticalParameterIndices() const
{
    auto combinedIndices = proportionalityFactorIndices;
    combinedIndices.insert(combinedIndices.end(), offsetParameterIndices.begin(), offsetParameterIndices.end());
    combinedIndices.insert(combinedIndices.end(), sigmaParameterIndices.begin(), sigmaParameterIndices.end());
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

void HierachicalOptimizationProblemWrapper::fillInitialParameters(gsl::span<double> buffer) const
{
    std::vector<double> full(wrappedProblem->costFun->numParameters());
    wrappedProblem->fillInitialParameters(full);
    fillFilteredParams(full, buffer);
}

void HierachicalOptimizationProblemWrapper::fillParametersMax(gsl::span<double> buffer) const
{
    std::vector<double> full(wrappedProblem->costFun->numParameters());
    wrappedProblem->fillParametersMax(full);
    fillFilteredParams(full, buffer);
}

void HierachicalOptimizationProblemWrapper::fillParametersMin(gsl::span<double> buffer) const
{
    std::vector<double> full(wrappedProblem->costFun->numParameters());
    wrappedProblem->fillParametersMin(full);
    fillFilteredParams(full, buffer);
}

void HierachicalOptimizationProblemWrapper::fillFilteredParams(const std::vector<double> &fullParams,
                                                               gsl::span<double> buffer) const
{
    auto hierarchical = dynamic_cast<HierachicalOptimizationWrapper *>(costFun.get());
    auto combinedIndices = hierarchical->getAnalyticalParameterIndices();
    parpe::fillFilteredParams(fullParams, combinedIndices, buffer);
}

std::unique_ptr<OptimizationReporter> HierachicalOptimizationProblemWrapper::getReporter() const {
    auto innerReporter = wrappedProblem->getReporter();
    auto outerReporter = std::make_unique<HierarchicalOptimizationReporter>(
                dynamic_cast<HierachicalOptimizationWrapper*>(costFun.get()),
                std::move(innerReporter->resultWriter));
    return outerReporter;
}

void fillFilteredParams(std::vector<double> const& valuesToFilter,
                        std::vector<int> const& sortedIndicesToExclude,
                        gsl::span<double> result)
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

double computeAnalyticalSigmas(int sigmaIdx, amici::AMICI_parameter_scaling scale,
                               const std::vector<std::vector<double> > &modelOutputsScaled,
                               const std::vector<std::vector<double> > &measurements,
                               AnalyticalParameterProvider &sigmaReader,
                               int numObservables, int numTimepoints)
{
    auto dependentConditions = sigmaReader.getConditionsForParameter(sigmaIdx);

    double enumerator = 0.0;
    double denominator = 0.0;

    for (auto const conditionIdx: dependentConditions) {
        auto dependentObservables = sigmaReader.getObservablesForParameter(sigmaIdx, conditionIdx);
        for(auto const observableIdx: dependentObservables) {
            for(int timeIdx = 0; timeIdx < numTimepoints; ++timeIdx) {
                double mes = measurements[conditionIdx][observableIdx + timeIdx * numObservables];
                if(!std::isnan(mes)) {
                    double scaledSim = modelOutputsScaled[conditionIdx][observableIdx + timeIdx * numObservables];
                    assert(!std::isnan(scaledSim));
                    enumerator += (mes - scaledSim) * (mes - scaledSim);
                    denominator += 1.0;
                }
            }
        }
    }

    double sigma = getScaledParameter(enumerator / denominator, scale);

    return sigma;
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



std::vector<double> spliceParameters(const gsl::span<double const> reducedParameters,
                                     const std::vector<int> &proportionalityFactorIndices,
                                     const std::vector<int> &offsetParameterIndices,
                                     const std::vector<int> &sigmaParameterIndices,
                                     const std::vector<double> &scalingFactors,
                                     const std::vector<double> &offsetParameters,
                                     const std::vector<double> &sigmaParameters) {

    std::vector<double> fullParameters(reducedParameters.size() + scalingFactors.size() + offsetParameters.size() + sigmaParameters.size());
    int idxScaling = 0;
    int idxOffset = 0;
    int idxSigma = 0;
    int idxRegular = 0;

    for(int i = 0; i < (signed) fullParameters.size(); ++i) {
        if((unsigned)idxScaling < proportionalityFactorIndices.size() && proportionalityFactorIndices[idxScaling] == i)
            fullParameters[i] = scalingFactors.at(idxScaling++);
        else if((unsigned)idxOffset < offsetParameterIndices.size() && offsetParameterIndices[idxOffset] == i)
            fullParameters[i] = offsetParameters.at(idxOffset++);
        else if((unsigned)idxSigma < sigmaParameterIndices.size() && sigmaParameterIndices[idxSigma] == i)
            fullParameters[i] = sigmaParameters.at(idxSigma++);
        else if((unsigned)idxRegular < reducedParameters.size())
            fullParameters[i] = reducedParameters.at(idxRegular++);
        else
            throw std::exception();
    }

    return fullParameters;
}


double computeNegLogLikelihood(std::vector <std::vector<double>> const& measurements,
                               const std::vector<std::vector<double> > &modelOutputsScaled,
                               std::vector <std::vector<double>> const& sigmas) {
    RELEASE_ASSERT(measurements.size() == modelOutputsScaled.size(), "");

    double nllh = 0.0;

    for (int conditionIdx = 0; (unsigned) conditionIdx < measurements.size(); ++conditionIdx) {
        nllh += computeNegLogLikelihood(measurements[conditionIdx], modelOutputsScaled[conditionIdx], sigmas[conditionIdx]);
    }

    return nllh;
}

double computeNegLogLikelihood(std::vector<double> const& measurements,
                               std::vector<double> const& modelOutputsScaled,
                               std::vector<double> const& sigmas) {
    double nllh = 0.0;
    constexpr double pi = atan(1)*4.0;

    RELEASE_ASSERT(measurements.size() == modelOutputsScaled.size(), "measurement/simulation output dimension mismatch");

    for(int i = 0; (unsigned) i < measurements.size(); ++i) {
        double mes = measurements[i];
        if(!std::isnan(mes)) {
            double sim = modelOutputsScaled[i];
            double sigmaSquared = sigmas[i] * sigmas[i];
            RELEASE_ASSERT(!std::isnan(sim), "");
            RELEASE_ASSERT(!std::isnan(sigmaSquared), "");
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

HierarchicalOptimizationReporter::HierarchicalOptimizationReporter(HierachicalOptimizationWrapper *gradFun, std::unique_ptr<OptimizationResultWriter> rw)
    : OptimizationReporter(gradFun, std::move(rw))
{
    hierarchicalWrapper = gradFun;
}

FunctionEvaluationStatus HierarchicalOptimizationReporter::evaluate(gsl::span<const double> parameters, double &fval, gsl::span<double> gradient) const
{
    if(beforeCostFunctionCall(parameters) != 0)
        return functionEvaluationFailure;

    if(gradient.data()) {
        if (!haveCachedGradient || !std::equal(parameters.begin(), parameters.end(),
                                               cachedParameters.begin())) {
            // Have to compute anew
            cachedStatus = hierarchicalWrapper->evaluate(parameters, cachedCost, cachedGradient,
                                                         cachedFullParameters, cachedFullGradient);
            haveCachedCost = true;
            haveCachedGradient = true;
        }
        // recycle old result
        std::copy(cachedGradient.begin(), cachedGradient.end(), gradient.begin());
        fval = cachedCost;
    } else {
        if (!haveCachedCost || !std::equal(parameters.begin(), parameters.end(),
                                           cachedParameters.begin())) {
            // Have to compute anew
            cachedStatus = hierarchicalWrapper->evaluate(parameters, cachedCost, gsl::span<double>(), cachedFullParameters, cachedFullGradient);
            haveCachedCost = true;
            haveCachedGradient = false;
        }
        fval = cachedCost;
    }

    // update cached parameters
    cachedParameters.resize(numParameters_);
    std::copy(parameters.begin(), parameters.end(), cachedParameters.begin());

    if(afterCostFunctionCall(parameters, cachedCost, gradient.data() ? cachedGradient : gsl::span<double>()) != 0)
        return functionEvaluationFailure;

    return cachedStatus;
}

void HierarchicalOptimizationReporter::finished(double optimalCost, gsl::span<const double> parameters, int exitStatus) const
{
    double timeElapsed = wallTimer.getTotal();

    if(cachedCost != optimalCost) {
        // the optimal value is not from the cached parameters and we did not get
        // the optimal full parameter vector. since we don't know them, rather set to nan
        cachedFullParameters.assign(cachedFullParameters.size(), NAN);
        std::copy(parameters.begin(), parameters.end(), cachedParameters.data());
    }

    cachedCost = optimalCost;

    if(resultWriter)
        resultWriter->saveLocalOptimizerResults(optimalCost, cachedFullParameters, timeElapsed, exitStatus);
}

bool HierarchicalOptimizationReporter::iterationFinished(gsl::span<const double> parameters, double objectiveFunctionValue, gsl::span<const double> objectiveFunctionGradient) const {
    double wallTimeIter = wallTimer.getRound(); //(double)(clock() - timeIterationBegin) / CLOCKS_PER_SEC;
    double wallTimeOptim = wallTimer.getTotal(); //double)(clock() - timeOptimizationBegin) / CLOCKS_PER_SEC;

    logmessage(LOGLVL_INFO, "iter: %d cost: %g time_iter: %gs time_optim: %gs", numIterations, objectiveFunctionValue, wallTimeIter, wallTimeOptim);

    if(resultWriter) {
        if(objectiveFunctionValue == cachedCost
                && (parameters.size() == 0 || std::equal(parameters.begin(), parameters.end(), cachedParameters.begin()))) {
            resultWriter->logLocalOptimizerIteration(numIterations, cachedFullParameters,
                                                     objectiveFunctionValue,
                                                     cachedFullGradient, // This might be misleading, the gradient could evaluated at other parameters if there was a line search inbetween
                                                     wallTimeIter);
        } else if (parameters.size()) {
            resultWriter->logLocalOptimizerIteration(numIterations, parameters,
                                                     objectiveFunctionValue,
                                                     objectiveFunctionGradient, // This might be misleading, the gradient could evaluated at other parameters if there was a line search inbetween
                                                     wallTimeIter);
        }
    }
    ++numIterations;

    return false;

}



} // namespace parpe
