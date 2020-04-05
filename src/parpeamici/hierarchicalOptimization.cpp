#include <parpeamici/hierarchicalOptimization.h>

#include <parpecommon/logging.h>
#include <parpeamici/multiConditionDataProvider.h>
#include <parpecommon/misc.h>
#include <parpecommon/parpeException.h>
#include <amici/misc.h>

#include <exception>
#include <cmath>

#ifndef __cpp_constexpr
// constexpr did not work on icc (ICC) 16.0.4 20160811
#define constexpr
#endif

namespace parpe {

HierarchicalOptimizationWrapper::HierarchicalOptimizationWrapper(
        std::unique_ptr<AmiciSummedGradientFunction> fun,
        int numConditions,
        int numObservables)
    : fun(std::move(fun)),
      numConditions(numConditions),
      numObservables(numObservables)
{
    scalingReader = std::make_unique<AnalyticalParameterHdf5Reader>();
    offsetReader  = std::make_unique<AnalyticalParameterHdf5Reader>();
    sigmaReader   = std::make_unique<AnalyticalParameterHdf5Reader>();

    if(fun)
        init();
}


HierarchicalOptimizationWrapper::HierarchicalOptimizationWrapper(
        std::unique_ptr<AmiciSummedGradientFunction > fun,
        H5::H5File const& file,
        std::string const& hdf5RootPath,
        int numConditions, int numObservables,
        ErrorModel errorModel)
    : fun(std::move(fun)),
      numConditions(numConditions),
      numObservables(numObservables),
      errorModel(errorModel)
{
    scalingReader = std::make_unique<AnalyticalParameterHdf5Reader>(
                file,
                hdf5RootPath + "/scalingParameterIndices",
                hdf5RootPath + "/scalingParametersMapToObservables");

    offsetReader = std::make_unique<AnalyticalParameterHdf5Reader>(
                file,
                hdf5RootPath + "/offsetParameterIndices",
                hdf5RootPath + "/offsetParametersMapToObservables");

    sigmaReader = std::make_unique<AnalyticalParameterHdf5Reader>(
                file,
                hdf5RootPath + "/sigmaParameterIndices",
                hdf5RootPath + "/sigmaParametersMapToObservables");

    init();
}


HierarchicalOptimizationWrapper::HierarchicalOptimizationWrapper(
        std::unique_ptr<AmiciSummedGradientFunction > fun,
        std::unique_ptr<AnalyticalParameterProvider> scalingReader,
        std::unique_ptr<AnalyticalParameterProvider> offsetReader,
        std::unique_ptr<AnalyticalParameterProvider> sigmaReader,
        int numConditions, int numObservables,
        ErrorModel errorModel)
    : fun(std::move(fun)),
      scalingReader(std::move(scalingReader)),
      offsetReader(std::move(offsetReader)),
      sigmaReader(std::move(sigmaReader)),
      numConditions(numConditions),
      numObservables(numObservables),
      errorModel(errorModel)
{
    init();
}


void HierarchicalOptimizationWrapper::init() {
    if(errorModel != ErrorModel::normal) {
        throw ParPEException("Only gaussian noise is supported so far.");
    }

    /* Some functions currently expect these lists to be sorted, therefore
     * ensure sorting right away (if sorting here, also need to reorder/reindex
     * scalingFactorIdx in mapping table -> difficult) */
    proportionalityFactorIndices =
            this->scalingReader->getOptimizationParameterIndices();
    RELEASE_ASSERT(
                std::is_sorted(this->proportionalityFactorIndices.begin(),
                               this->proportionalityFactorIndices.end()), "");

    offsetParameterIndices =
            this->offsetReader->getOptimizationParameterIndices();
    RELEASE_ASSERT(
                std::is_sorted(this->offsetParameterIndices.begin(),
                               this->offsetParameterIndices.end()), "");

    sigmaParameterIndices =
            this->sigmaReader->getOptimizationParameterIndices();
    RELEASE_ASSERT(std::is_sorted(this->sigmaParameterIndices.begin(),
                                  this->sigmaParameterIndices.end()), "");

    if(fun) {
        std::stringstream ss;
        ss<<"HierarchicalOptimizationWrapper parameters: "
           <<fun->numParameters()<<" total, "
           <<numParameters()<< " numerical, "
           <<proportionalityFactorIndices.size()<<" proportionality, "
           <<offsetParameterIndices.size()<<" offset, "
           <<sigmaParameterIndices.size()<<" sigma\n";
        Logger logger;
        logger.logmessage(LOGLVL_DEBUG, ss.str());
    }
}


FunctionEvaluationStatus HierarchicalOptimizationWrapper::evaluate(
        gsl::span<const double> parameters,
        double &fval,
        gsl::span<double> gradient,
        Logger *logger,
        double *cpuTime) const {

    std::vector<double> fullParameters;
    std::vector<double> fullGradient;

    return evaluate(parameters, fval, gradient, fullParameters,
                    fullGradient, logger, cpuTime);
}

FunctionEvaluationStatus
HierarchicalOptimizationWrapper::evaluate(
        gsl::span<const double> reducedParameters,
        double &fval,
        gsl::span<double> gradient,
        std::vector<double> &fullParameters,
        std::vector<double> &fullGradient, Logger *logger, double *cpuTime) const
{
    WallTimer walltimer;
    FunctionEvaluationStatus status;

    if(reducedParameters.size() != (unsigned)numParameters()) {
        throw ParPEException("Reduced parameter vector size "
                             + std::to_string(reducedParameters.size())
                             + " does not match numParameters "
                             + std::to_string(numParameters()));
    }
    RELEASE_ASSERT(gradient.empty()
                   || gradient.size() == reducedParameters.size(), "");

    if(numProportionalityFactors() == 0
            && numOffsetParameters() == 0
            && numSigmaParameters() == 0) {
        // nothing to do, just pass through

        // evaluate for all conditions
        std::vector<int> dataIndices(numConditions);
        std::iota(dataIndices.begin(), dataIndices.end(), 0);
        return fun->evaluate(reducedParameters, dataIndices,
                             fval, gradient, logger, cpuTime);
    }


    // evaluate with scaling parameters set to 1 and offsets to 0
    std::vector<std::vector<double> > modelOutput;
    try {
        modelOutput = getUnscaledModelOutputs(reducedParameters, logger, cpuTime);
    } catch (ParPEException const &e) {
        return FunctionEvaluationStatus::functionEvaluationFailure;
    }

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

    if(logger) {
        std::stringstream ss;
        ss<<"scalings "<<scalings;
        logger->logmessage(LOGLVL_DEBUG, ss.str());
        ss.str(std::string());
        ss<<"sigmas "<<sigmas;
        logger->logmessage(LOGLVL_DEBUG, ss.str());
    }
    // splice parameter vector we get from optimizer with analytically
    // computed parameters
    fullParameters = spliceParameters(
                reducedParameters, proportionalityFactorIndices,
                offsetParameterIndices, sigmaParameterIndices,
                scalings, offsets, sigmas);


    // evaluate with analytical scaling parameters
    double cpuTimeInner = 0.0;
    status = evaluateWithOptimalParameters(fullParameters, sigmas,
                                           measurements, modelOutput,
                                           fval, gradient,
                                           fullGradient, logger, &cpuTimeInner);

    if(cpuTime)
        *cpuTime += cpuTimeInner + walltimer.getTotal();

    return status;
}


std::vector<double>
HierarchicalOptimizationWrapper::getDefaultScalingFactors() const
{
    auto result = std::vector<double>(numProportionalityFactors());

    for(int i = 0; i < numProportionalityFactors(); ++i) {
        result[i] = getDefaultScalingFactor(
                    fun->getParameterScaling(proportionalityFactorIndices[i]));
    }

    return result;
}


std::vector<double>
HierarchicalOptimizationWrapper::getDefaultOffsetParameters() const
{
    auto result = std::vector<double>(numOffsetParameters());

    for(int i = 0; i < numOffsetParameters(); ++i) {
        result[i] = getDefaultOffsetParameter(
                    fun->getParameterScaling(offsetParameterIndices[i]));
    }

    return result;
}


std::vector<double>
HierarchicalOptimizationWrapper::getDefaultSigmaParameters() const
{
    auto result = std::vector<double>(numSigmaParameters());

    for(int i = 0; i < numSigmaParameters(); ++i) {
        // default sigma is the same as default scaling
        result[i] = getDefaultScalingFactor(
                    fun->getParameterScaling(sigmaParameterIndices[i]));
    }

    return result;

}

std::vector<std::vector<double> >
HierarchicalOptimizationWrapper::getUnscaledModelOutputs(
        const gsl::span<const double> reducedParameters, Logger *logger,
        double *cpuTime) const
{
    // run simulations, collect outputs
    auto scalingDummy = getDefaultScalingFactors();
    auto offsetDummy = getDefaultOffsetParameters();
    auto sigmaDummy = getDefaultSigmaParameters();

    // splice hidden scaling parameter and external parameters
    auto fullParameters = spliceParameters(
                reducedParameters, proportionalityFactorIndices,
                offsetParameterIndices, sigmaParameterIndices,
                scalingDummy, offsetDummy, sigmaDummy);

    std::vector<std::vector<double> > modelOutput(numConditions);
    auto status = fun->getModelOutputs(fullParameters, modelOutput,
                                       logger, cpuTime);
    if(status != FunctionEvaluationStatus::functionEvaluationSuccess)
        throw ParPEException("Function evaluation failed.");

    return modelOutput;
}

std::vector<double> HierarchicalOptimizationWrapper::computeAnalyticalScalings(
        std::vector<std::vector<double>> const& measurements,
        std::vector<std::vector<double>> const& modelOutputsUnscaled) const
{
    int numProportionalityFactors = proportionalityFactorIndices.size();
    std::vector<double> proportionalityFactors(numProportionalityFactors);

    for(int scalingIdx = 0; scalingIdx < numProportionalityFactors;
        ++scalingIdx) {

        auto proportionalityFactor = parpe::computeAnalyticalScalings(
                    scalingIdx, modelOutputsUnscaled, measurements,
                    *scalingReader, numObservables);
        auto scale = fun->getParameterScaling(
                    proportionalityFactorIndices[scalingIdx]);
        proportionalityFactors[scalingIdx] =
                getScaledParameter(proportionalityFactor, scale);
    }

    return proportionalityFactors;
}

void HierarchicalOptimizationWrapper::applyOptimalScalings(
        std::vector<double> const& proportionalityFactors,
        std::vector<std::vector<double> > &modelOutputs) const {

    for(int i = 0; (unsigned) i < proportionalityFactors.size(); ++i) {
        double scaling = getUnscaledParameter(
                    proportionalityFactors[i],
                    fun->getParameterScaling(proportionalityFactorIndices[i]));

        applyOptimalScaling(i, scaling, modelOutputs,
                            *scalingReader, numObservables);
    }
}


std::vector<double> HierarchicalOptimizationWrapper::computeAnalyticalOffsets(
        std::vector<std::vector<double>> const& measurements,
        std::vector<std::vector<double> > &modelOutputsUnscaled) const
{
    int numOffsetParameters = offsetParameterIndices.size();
    std::vector<double> offsetParameters(numOffsetParameters);

    for(int i = 0; i < numOffsetParameters; ++i) {
        auto offsetParameter = parpe::computeAnalyticalOffsets(
                    i, modelOutputsUnscaled, measurements,
                    *offsetReader, numObservables);
        auto scale = fun->getParameterScaling(offsetParameterIndices[i]);
        offsetParameters[i] = getScaledParameter(offsetParameter, scale);
    }

    return offsetParameters;
}

std::vector<double> HierarchicalOptimizationWrapper::computeAnalyticalSigmas(
        const std::vector<std::vector<double> > &measurements,
        std::vector<std::vector<double> > const &modelOutputsScaled) const
{
    int numSigmas = sigmaParameterIndices.size();
    std::vector<double> sigmas(numSigmas);

    for(int i = 0; i < numSigmas; ++i) {
        auto sigma = parpe::computeAnalyticalSigmas(
                    i,
                    modelOutputsScaled, measurements,
                    *sigmaReader, numObservables);
        auto scale = fun->getParameterScaling(sigmaParameterIndices[i]);
        sigmas[i] = getScaledParameter(sigma, scale);
    }
    return sigmas;
}

void HierarchicalOptimizationWrapper::applyOptimalOffsets(
        std::vector<double> const& offsetParameters,
        std::vector<std::vector<double> > &modelOutputs) const {

    for(int i = 0; (unsigned) i < offsetParameters.size(); ++i) {
        double offset = getUnscaledParameter(
                    offsetParameters[i],
                    fun->getParameterScaling(offsetParameterIndices[i]));
        applyOptimalOffset(i, offset, modelOutputs,
                           *offsetReader, numObservables);
    }
}

void HierarchicalOptimizationWrapper::fillInAnalyticalSigmas(
        std::vector<std::vector<double> > &allSigmas,
        std::vector<double> const& analyticalSigmas) const
{
    for(int sigmaParameterIdx = 0;
        (unsigned) sigmaParameterIdx < analyticalSigmas.size();
        ++sigmaParameterIdx) {

        // sigma value will be used for likelihood computation
        // and not passed to AMICI -> unscale
        auto sigmaParameterValue = getUnscaledParameter(
                    analyticalSigmas[sigmaParameterIdx],
                    fun->getParameterScaling(
                        sigmaParameterIndices[sigmaParameterIdx]));

        auto dependentConditions =
                sigmaReader->getConditionsForParameter(sigmaParameterIdx);

        for (auto const conditionIdx: dependentConditions) {

            int numTimepoints = allSigmas[conditionIdx].size() / numObservables;

            auto dependentObservables =
                    sigmaReader->getObservablesForParameter(
                        sigmaParameterIdx, conditionIdx);

            for(auto const observableIdx: dependentObservables) {

                RELEASE_ASSERT(observableIdx < numObservables, "");

                for(int timeIdx = 0; timeIdx < numTimepoints; ++timeIdx) {
                    // NOTE: this must be in sync with data ordering in AMICI
                    // (assumes row-major)
                    RELEASE_ASSERT(
                                std::isnan(
                                    allSigmas[conditionIdx][observableIdx + timeIdx * numObservables]),
                            "Expected NaN value for sigma parameters being "
                            "estimated, but got non-NAN.");
                    allSigmas[conditionIdx][observableIdx + timeIdx * numObservables] = sigmaParameterValue;
                }
            }
        }
    }
}


FunctionEvaluationStatus
HierarchicalOptimizationWrapper::evaluateWithOptimalParameters(
        std::vector<double> const& fullParameters,
        std::vector<double> const& sigmas,
        std::vector<std::vector<double>> const& measurements,
        std::vector<std::vector<double>> const& modelOutputsScaled,
        double &fval,
        const gsl::span<double> gradient,
        std::vector<double>& fullGradient,
        Logger *logger, double *cpuTime) const {

    if(!gradient.empty()) {
        fval = NAN;
        // simulate with updated theta for sensitivities
        // simulate all datasets
        std::vector<int> dataIndices(numConditions);
        std::iota(dataIndices.begin(), dataIndices.end(), 0);

        // Need intermediary buffer because optimizer expects
        // fewer parameters than `fun` delivers
        fullGradient.resize(fullParameters.size());
        auto status = fun->evaluate(fullParameters, dataIndices, fval,
                                    fullGradient, logger, cpuTime);
        if(status != functionEvaluationSuccess)
            return status;

        // Filter gradient for those parameters expected by the optimizer
        auto analyticalParameterIndices = getAnalyticalParameterIndices();
        fillFilteredParams(fullGradient, analyticalParameterIndices, gradient);

        // Check if gradient w.r.t. analytical parameters is 0
        checkGradientForAnalyticalParameters(
                    fullGradient, analyticalParameterIndices, 1e-8);

    } else {
        auto fullSigmaMatrices = fun->getAllSigmas();
        if(!sigmaParameterIndices.empty()) {
            fillInAnalyticalSigmas(fullSigmaMatrices, sigmas);
        }

        // ... to compute negative log-likelihood
        fval = computeNegLogLikelihood(measurements, modelOutputsScaled,
                                       fullSigmaMatrices);
    }

    return std::isfinite(fval) ?
                functionEvaluationSuccess : functionEvaluationFailure;
}


int HierarchicalOptimizationWrapper::numParameters() const {
    return fun->numParameters() - numProportionalityFactors()
            - numOffsetParameters() - numSigmaParameters();
}

int HierarchicalOptimizationWrapper::numProportionalityFactors() const {
    return proportionalityFactorIndices.size();
}

const std::vector<int> &
HierarchicalOptimizationWrapper::getProportionalityFactorIndices() const
{
    return proportionalityFactorIndices;
}


AnalyticalParameterHdf5Reader::AnalyticalParameterHdf5Reader(
        H5::H5File const& file,
        std::string analyticalParameterIndicesPath,
        std::string mapPath)
    : mapPath(std::move(mapPath)),
      analyticalParameterIndicesPath(std::move(analyticalParameterIndicesPath))
{
    auto lock = hdf5MutexGetLock();
    this->file = file; // copy while mutex is locked!
    readParameterConditionObservableMappingFromFile();
}


int HierarchicalOptimizationWrapper::numOffsetParameters() const {
    return offsetParameterIndices.size();
}

int HierarchicalOptimizationWrapper::numSigmaParameters() const
{
    return sigmaParameterIndices.size();
}

const std::vector<int> &
HierarchicalOptimizationWrapper::getOffsetParameterIndices() const
{
    return offsetParameterIndices;
}

const std::vector<int> &
HierarchicalOptimizationWrapper::getSigmaParameterIndices() const
{
    return sigmaParameterIndices;
}


std::vector<int>
HierarchicalOptimizationWrapper::getAnalyticalParameterIndices() const
{
    auto combinedIndices = proportionalityFactorIndices;
    combinedIndices.insert(combinedIndices.end(),
                           offsetParameterIndices.begin(),
                           offsetParameterIndices.end());
    combinedIndices.insert(combinedIndices.end(),
                           sigmaParameterIndices.begin(),
                           sigmaParameterIndices.end());
    std::sort(combinedIndices.begin(), combinedIndices.end());

    return combinedIndices;
}


std::vector<int>
AnalyticalParameterHdf5Reader::getConditionsForParameter(
        int parameterIndex) const {
    std::vector<int> result;
    result.reserve(mapping[parameterIndex].size());
    for (auto const& kvp : mapping[parameterIndex])
        result.push_back(kvp.first);
    return result;
}


const std::vector<int> &
AnalyticalParameterHdf5Reader::getObservablesForParameter(
        int parameterIndex, int conditionIdx) const {
    return mapping[parameterIndex].at(conditionIdx);
}


std::vector<int>
AnalyticalParameterHdf5Reader::getOptimizationParameterIndices() const {
    auto lock = hdf5MutexGetLock();
    std::vector<int> analyticalParameterIndices;

    if(file.nameExists(analyticalParameterIndicesPath)) {
        auto dataset = file.openDataSet(analyticalParameterIndicesPath);
        auto dataspace = dataset.getSpace();

        auto ndims = dataspace.getSimpleExtentNdims();
        if(ndims != 1)
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

int AnalyticalParameterHdf5Reader::getNumAnalyticalParameters() const
{
    hsize_t numAnalyticalParameters = 0;
    auto lock = hdf5MutexGetLock();

    if(file.nameExists(analyticalParameterIndicesPath)) {
        auto dataset = file.openDataSet(analyticalParameterIndicesPath);
        auto dataspace = dataset.getSpace();
        auto ndims = dataspace.getSimpleExtentNdims();
        if(ndims != 1)
            throw ParPEException(
                    "Invalid dimension in getOptimizationParameterIndices.");
        dataspace.getSimpleExtentDims(&numAnalyticalParameters);
    }

    return numAnalyticalParameters;
}

void
AnalyticalParameterHdf5Reader::readParameterConditionObservableMappingFromFile()
{
    auto lock = hdf5MutexGetLock();
    H5_SAVE_ERROR_HANDLER;
    try {
        int numScalings = getNumAnalyticalParameters();
        auto dataset = file.openDataSet(mapPath);
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
    } catch (H5::FileIException&) {
        return;
    }
    H5_RESTORE_ERROR_HANDLER;

}

std::vector<int> AnalyticalParameterHdf5Reader::readRawMap(
        H5::DataSet& dataset, hsize_t& nRows, hsize_t& nCols)
{
    auto dataspace = dataset.getSpace();
    auto ndims = dataspace.getSimpleExtentNdims();
    if(ndims != 2)
        throw ParPEException(
                "Invalid dimension for analytical parameter map, expected 2.");

    hsize_t dims[ndims];
    dataspace.getSimpleExtentDims(dims);
    nRows = dims[0];
    nCols = dims[1];
    if(nRows && nCols != 3)
        throw ParPEException(
                "Invalid dimension for analytical parameter map, "
                "expected 3 columns.");

    std::vector<int> rawMap(nRows * nCols);
    dataset.read(rawMap.data(), H5::PredType::NATIVE_INT);

    return rawMap;
}

HierarchicalOptimizationProblemWrapper::HierarchicalOptimizationProblemWrapper(
        std::unique_ptr<OptimizationProblem> problemToWrap,
        const MultiConditionDataProviderHDF5 *dataProvider)
    : wrapped_problem_(std::move(problemToWrap))
{
    logger_ = std::make_unique<Logger>(*wrapped_problem_->logger_);
    auto wrappedFun =
            dynamic_cast<SummedGradientFunctionGradientFunctionAdapter<int>*>(
                wrapped_problem_->cost_fun_.get());

    auto model = dataProvider->getModel();

    auto lock = hdf5MutexGetLock();
    cost_fun_.reset(
                new HierarchicalOptimizationWrapper(
                      std::unique_ptr<AmiciSummedGradientFunction>(
                          dynamic_cast<AmiciSummedGradientFunction*>(
                            wrappedFun->getWrappedFunction())),
                      dataProvider->getHdf5FileId(), "/",
                      dataProvider->getNumberOfSimulationConditions(),
                      model->nytrue,
                      ErrorModel::normal));
}

HierarchicalOptimizationProblemWrapper::HierarchicalOptimizationProblemWrapper(
        std::unique_ptr<OptimizationProblem> problemToWrap,
        std::unique_ptr<HierarchicalOptimizationWrapper> costFun,
        std::unique_ptr<Logger> logger)
    : OptimizationProblem(std::move(costFun),
                          std::move(logger)),
      wrapped_problem_(std::move(problemToWrap))
{

}

HierarchicalOptimizationProblemWrapper
::~HierarchicalOptimizationProblemWrapper()
{
    // Avoid double delete.
    // This will be destroyed when wrappedProblem goes out of scope!
    dynamic_cast<HierarchicalOptimizationWrapper *>(
                cost_fun_.get())->fun.release();
}

void HierarchicalOptimizationProblemWrapper::fillInitialParameters(
        gsl::span<double> buffer) const
{
    std::vector<double> full(wrapped_problem_->cost_fun_->numParameters());
    wrapped_problem_->fillInitialParameters(full);
    fillFilteredParams(full, buffer);
}

void HierarchicalOptimizationProblemWrapper::fillParametersMax(
        gsl::span<double> buffer) const
{
    std::vector<double> full(wrapped_problem_->cost_fun_->numParameters());
    wrapped_problem_->fillParametersMax(full);
    fillFilteredParams(full, buffer);
}

void HierarchicalOptimizationProblemWrapper::fillParametersMin(
        gsl::span<double> buffer) const
{
    std::vector<double> full(wrapped_problem_->cost_fun_->numParameters());
    wrapped_problem_->fillParametersMin(full);
    fillFilteredParams(full, buffer);
}

void HierarchicalOptimizationProblemWrapper::fillFilteredParams(
        const std::vector<double> &fullParams,
        gsl::span<double> buffer) const
{
    auto hierarchical = dynamic_cast<HierarchicalOptimizationWrapper *>(
                cost_fun_.get());
    auto combinedIndices = hierarchical->getAnalyticalParameterIndices();
    parpe::fillFilteredParams(fullParams, combinedIndices, buffer);
}

std::unique_ptr<OptimizationReporter>
HierarchicalOptimizationProblemWrapper::getReporter() const {
    auto innerReporter = wrapped_problem_->getReporter();
    auto outerReporter = std::unique_ptr<OptimizationReporter>(
                new HierarchicalOptimizationReporter(
                    dynamic_cast<HierarchicalOptimizationWrapper*>(cost_fun_.get()),
                    std::move(innerReporter->result_writer_),
                    std::make_unique<Logger>(*logger_)
                    ));
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
    RELEASE_ASSERT(resultIdx ==
                   (unsigned) valuesToFilter.size()
                   - sortedIndicesToExclude.size(),
                   "");
}

double getDefaultScalingFactor(amici::ParameterScaling scaling)
{
    switch (scaling) {
    case amici::ParameterScaling::none:
        return 1.0;
    case amici::ParameterScaling::log10:
        return 0.0;
    default:
        throw ParPEException(
                    "Parameter scaling must be ParameterScaling::log10 "
                    "or ParameterScaling::none.");
    }
}

double getDefaultOffsetParameter(amici::ParameterScaling scaling)
{
    switch (scaling) {
    case amici::ParameterScaling::none:
        return 0.0;
    case amici::ParameterScaling::log10:
        return -std::numeric_limits<double>::infinity();
    default:
        throw ParPEException(
                    "Parameter scaling must be ParameterScaling::log10 "
                    "or ParameterScaling::none.");
    }
}


double computeAnalyticalScalings(
        int scalingIdx,
        const std::vector<std::vector<double> > &modelOutputsUnscaled,
        const std::vector<std::vector<double> > &measurements,
        AnalyticalParameterProvider const& scalingReader,
        int numObservables) {

    auto dependentConditions =
            scalingReader.getConditionsForParameter(scalingIdx);

    double enumerator = 0.0;
    double denominator = 0.0;

    for (auto const conditionIdx: dependentConditions) {
        auto dependentObservables =
                scalingReader.getObservablesForParameter(
                    scalingIdx, conditionIdx);
        int numTimepoints = measurements[conditionIdx].size() / numObservables;

        for(auto const observableIdx: dependentObservables) {
            for(int timeIdx = 0; timeIdx < numTimepoints; ++timeIdx) {

                double mes = measurements[conditionIdx][observableIdx + timeIdx * numObservables];
                if(!std::isnan(mes)) {
                    // NOTE: this must be in sync with data ordering in AMICI (assumes row-major)
                    double sim = modelOutputsUnscaled[conditionIdx][observableIdx + timeIdx * numObservables];
                    // std::cout<<scalingIdx<<"\t"<<conditionIdx<<"\t"<<observableIdx<<"\t"<<timeIdx<<"\t"<<mes<<"\t"<<sim<<std::endl;
                    if(std::isnan(sim)) {
                        logmessage(LOGLVL_WARNING,
                                   "In computeAnalyticalScalings %d: "
                                   "Simulation is NaN for condition %d "
                                   "observable %d timepoint %d", scalingIdx,
                                   conditionIdx, observableIdx, timeIdx);
                    }
                    if(sim < 0 && sim > -1e-18) {
                        // negative values due to numerical errors
                        // TODO: some outputs may be validly < 0
                        logmessage(LOGLVL_WARNING,
                                   "In computeAnalyticalScalings %d: "
                                   "Simulation is %g < 0 for condition %d "
                                   "observable %d timepoint %d. "
                                   "Setting to 0.0.", scalingIdx, sim,
                                   conditionIdx, observableIdx, timeIdx);
                        sim = 0.0;
                    }

                    enumerator += sim * mes;
                    denominator += sim * sim;
                }
            }
        }
    }

    if(denominator == 0.0) {
        logmessage(LOGLVL_WARNING,
                   "In computeAnalyticalScalings: denominator is 0.0 for "
                   "scaling parameter " + std::to_string(scalingIdx)
                   + ". Probably model output is always 0.0 and scaling, "
                     "thus, not used. Setting scaling parameter to 1.0.");
        return 1.0;
    }
    double scaling = enumerator / denominator;

    return scaling;
//    constexpr double upper_bound = 1e10;

//    // too large values of scaling parameters cause problems in backwards
//    // integration for adjoint sensitivities
//    if(upper_bound > scaling)
//        return scaling;

//    logmessage(LOGLVL_WARNING,
//               "In computeAnalyticalScalings: force-bounding scaling parameter "
//               + std::to_string(scalingIdx) + " which was "
//               + std::to_string(scaling) + " to "
//               + std::to_string(upper_bound));
//    return upper_bound;
}


double computeAnalyticalOffsets(
        int offsetIdx,
        std::vector<std::vector<double>> const& modelOutputsUnscaled,
        std::vector<std::vector<double>> const& measurements,
        AnalyticalParameterProvider& offsetReader,
        int numObservables)
{
    auto dependentConditions =
            offsetReader.getConditionsForParameter(offsetIdx);

    double enumerator = 0.0;
    double denominator = 0.0;

    for (auto const conditionIdx: dependentConditions) {
        auto dependentObservables =
                offsetReader.getObservablesForParameter(
                    offsetIdx, conditionIdx);
        int numTimepoints = measurements[conditionIdx].size() / numObservables;
        for(auto const observableIdx: dependentObservables) {
            for(int timeIdx = 0; timeIdx < numTimepoints; ++timeIdx) {
                double mes =
                        measurements[conditionIdx][observableIdx + timeIdx * numObservables];
                if(!std::isnan(mes)) {
                    double sim =
                            modelOutputsUnscaled[conditionIdx][observableIdx + timeIdx * numObservables];
                    if(std::isnan(sim)) {
                        logmessage(LOGLVL_WARNING,
                                   "In computeAnalyticalOffsets %d: "
                                   "Simulation is NaN for condition %d "
                                   "observable %d timepoint %d", offsetIdx,
                                   conditionIdx, observableIdx, timeIdx);
                    }
                    enumerator += mes - sim;
                    denominator += 1.0;
                }
            }
        }
    }

    if(denominator == 0.0) {
        logmessage(LOGLVL_WARNING,
                   "In computeAnalyticalOffsets: denominator is 0.0 "
                   "for offset parameter " + std::to_string(offsetIdx)
                   + ". This probably means that there exists no measurement "
                     "using this parameter. Setting offset to 0.0.");
        return 0.0;
    }

    return enumerator / denominator;
}

double computeAnalyticalSigmas(
        int sigmaIdx,
        const std::vector<std::vector<double> > &modelOutputsScaled,
        const std::vector<std::vector<double> > &measurements,
        AnalyticalParameterProvider const& sigmaReader,
        int numObservables,
        double epsilonAbs, double epsilonRel)
{
    auto dependentConditions = sigmaReader.getConditionsForParameter(sigmaIdx);

    double enumerator = 0.0;
    double denominator = 0.0;

    double maxAbsMeasurement = 0.0;

    for (auto const conditionIdx: dependentConditions) {
        auto dependentObservables =
                sigmaReader.getObservablesForParameter(sigmaIdx, conditionIdx);
        int numTimepoints = measurements[conditionIdx].size() / numObservables;
        for(auto const observableIdx: dependentObservables) {
            if(observableIdx >= numObservables) {
                throw ParPEException("computeAnalyticalSigmas: Invalid "
                                     "observableIdx >= numObservables.");
            }

            for(int timeIdx = 0; timeIdx < numTimepoints; ++timeIdx) {
                double mes = measurements[conditionIdx][observableIdx + timeIdx * numObservables];
                if(!std::isnan(mes)) {
                    double scaledSim = modelOutputsScaled[conditionIdx][observableIdx + timeIdx * numObservables];
                    // std::cout<<sigmaIdx<<"\t"<<conditionIdx<<"\t"<<observableIdx<<"\t"<<timeIdx<<"\t"<<mes<<"\t"<<scaledSim<<std::endl;
                    if(std::isnan(scaledSim)) {
                        logmessage(LOGLVL_WARNING,
                                   "In computeAnalyticalSigmas %d: "
                                   "Simulation is NaN for condition %d observable %d timepoint %d",
                                   sigmaIdx, conditionIdx, observableIdx, timeIdx);
                    }
                    enumerator += (mes - scaledSim) * (mes - scaledSim);
                    denominator += 1.0;

                    maxAbsMeasurement = std::max(maxAbsMeasurement, std::abs(mes));
                }
            }
        }
    }

    if(denominator == 0.0) {
        logmessage(LOGLVL_WARNING,
                   "In computeAnalyticalSigmas: Denominator is 0.0 for sigma parameter "
                   + std::to_string(sigmaIdx)
                   + ". This probably means that there exists no measurement using this parameter.");
    }

    double sigma = std::sqrt(enumerator / denominator);

    epsilonAbs = std::max(epsilonRel * maxAbsMeasurement, epsilonAbs);

    if(sigma < epsilonAbs) {
        // Must not return sigma = 0.0
        logmessage(LOGLVL_WARNING, "In computeAnalyticalSigmas " + std::to_string(sigmaIdx)
                   + ": Computed sigma < epsilon. Setting to " + std::to_string(epsilonAbs));
        return epsilonAbs;
    }

    return sigma;
}


void applyOptimalScaling(int scalingIdx, double scalingLin,
                         std::vector<std::vector<double> > &modelOutputs,
                         AnalyticalParameterProvider const& scalingReader,
                         int numObservables) {
    auto dependentConditions = scalingReader.getConditionsForParameter(scalingIdx);
    for (auto const conditionIdx: dependentConditions) {
        int numTimepoints = modelOutputs[conditionIdx].size() / numObservables;
        auto dependentObservables = scalingReader.getObservablesForParameter(scalingIdx, conditionIdx);
        for(auto const observableIdx: dependentObservables) {
            if(observableIdx >= numObservables) {
                throw ParPEException("applyOptimalOffset: Invalid observableIdx >= numObservables.");
            }

            for(int timeIdx = 0; timeIdx < numTimepoints; ++timeIdx) {
                // NOTE: this must be in sync with data ordering in AMICI (assumes row-major)
                modelOutputs[conditionIdx][observableIdx + timeIdx * numObservables] *= scalingLin;
            }
        }
    }
}



void applyOptimalOffset(int offsetIdx, double offsetLin,
                        std::vector<std::vector<double> > &modelOutputs,
                        const AnalyticalParameterProvider &offsetReader,
                        int numObservables) {
    auto dependentConditions = offsetReader.getConditionsForParameter(offsetIdx);
    for (auto const conditionIdx: dependentConditions) {
        int numTimepoints = modelOutputs[conditionIdx].size() / numObservables;
        auto dependentObservables = offsetReader.getObservablesForParameter(offsetIdx, conditionIdx);
        for(auto const observableIdx: dependentObservables) {
            if(observableIdx >= numObservables) {
                throw ParPEException("applyOptimalOffset: Invalid observableIdx >= numObservables.");
            }
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

    std::vector<double> fullParameters(
                reducedParameters.size() + scalingFactors.size()
                + offsetParameters.size() + sigmaParameters.size());
    int idxScaling = 0;
    int idxOffset = 0;
    int idxSigma = 0;
    int idxRegular = 0;

    for(int i = 0; i < (signed) fullParameters.size(); ++i) {
        if((unsigned)idxScaling < proportionalityFactorIndices.size()
                && proportionalityFactorIndices[idxScaling] == i)
            fullParameters[i] = scalingFactors.at(idxScaling++);
        else if((unsigned)idxOffset < offsetParameterIndices.size()
                && offsetParameterIndices[idxOffset] == i)
            fullParameters[i] = offsetParameters.at(idxOffset++);
        else if((unsigned)idxSigma < sigmaParameterIndices.size()
                && sigmaParameterIndices[idxSigma] == i)
            fullParameters[i] = sigmaParameters.at(idxSigma++);
        else if((unsigned)idxRegular < reducedParameters.size())
            fullParameters[i] = reducedParameters.at(idxRegular++);
        else
            throw std::exception();
    }

    RELEASE_ASSERT((unsigned) idxScaling == proportionalityFactorIndices.size(),
                   "")
    RELEASE_ASSERT((unsigned) idxOffset == offsetParameterIndices.size(), "")
    RELEASE_ASSERT((unsigned) idxSigma == sigmaParameterIndices.size(), "")
    RELEASE_ASSERT((unsigned) idxRegular == reducedParameters.size(), "")

    return fullParameters;
}


double computeNegLogLikelihood(
        std::vector<std::vector<double>> const& measurements,
        std::vector<std::vector<double>> const& modelOutputsScaled,
        std::vector<std::vector<double>> const& sigmas) {
    RELEASE_ASSERT(measurements.size() == modelOutputsScaled.size(), "");

    double nllh = 0.0;

    for (int conditionIdx = 0;
         (unsigned) conditionIdx < measurements.size(); ++conditionIdx) {
        nllh += computeNegLogLikelihood
                (measurements[conditionIdx], modelOutputsScaled[conditionIdx],
                 sigmas[conditionIdx]);
        if(std::isnan(nllh))
            return nllh;
    }

    return nllh;
}

double computeNegLogLikelihood(std::vector<double> const& measurements,
                               std::vector<double> const& modelOutputsScaled,
                               std::vector<double> const& sigmas) {
    double nllh = 0.0;

    RELEASE_ASSERT(measurements.size() == modelOutputsScaled.size(),
                   "measurement/simulation output dimension mismatch");

    for(int i = 0; (unsigned) i < measurements.size(); ++i) {
        double mes = measurements[i];
        if(!std::isnan(mes)) {
            double sim = modelOutputsScaled[i];
            double sigmaSquared = sigmas[i] * sigmas[i];
            if(std::isnan(sim)) {
                logmessage(LOGLVL_WARNING, "Simulation is NaN for data point %d", i);
                return std::numeric_limits<double>::quiet_NaN();
            }
            if(std::isnan(sigmaSquared)) {
                logmessage(LOGLVL_WARNING, "Sigma is NaN for data point %d", i);
                return std::numeric_limits<double>::quiet_NaN();
            }
            if(sigmaSquared < 0.0) {
                logmessage(LOGLVL_WARNING, "Negative sigma for data point %d", i);
                return std::numeric_limits<double>::quiet_NaN();
            }

            double diff = mes - sim;
            diff *= diff;
            nllh += log(2.0 * M_PI * sigmaSquared) + diff / sigmaSquared;
        }
    }

    nllh /= 2.0;
    return nllh;
}

std::vector<int>
AnalyticalParameterProviderDefault::getConditionsForParameter(
        int parameterIndex) const {
    return conditionsForParameter[parameterIndex];
}

const std::vector<int> &
AnalyticalParameterProviderDefault::getObservablesForParameter(
        int parameterIndex, int conditionIdx) const {
    return mapping[parameterIndex].at(conditionIdx);
}

std::vector<int> AnalyticalParameterProviderDefault::getOptimizationParameterIndices() const {
    return optimizationParameterIndices;
}

HierarchicalOptimizationReporter::HierarchicalOptimizationReporter(
        HierarchicalOptimizationWrapper *gradFun,
        std::unique_ptr<OptimizationResultWriter> rw,
        std::unique_ptr<Logger> logger)
    : OptimizationReporter(gradFun, std::move(rw), std::move(logger))
{
    hierarchical_wrapper_ = gradFun;
}

FunctionEvaluationStatus HierarchicalOptimizationReporter::evaluate(
        gsl::span<const double> parameters,
        double &fval, gsl::span<double> gradient, Logger *logger,
        double *cpuTime) const
{
    double myCpuTimeSec = 0.0;
    if(cpuTime)
        *cpuTime = 0.0;

    if(beforeCostFunctionCall(parameters) != 0)
        return functionEvaluationFailure;

    if(gradient.data()) {
        if (!have_cached_gradient_ || !std::equal(parameters.begin(), parameters.end(),
                                               cached_parameters_.begin())) {
            // Have to compute anew
            cached_status_ = hierarchical_wrapper_->evaluate(
                        parameters, cached_cost_, cached_gradient_,
                        cached_full_parameters_, cached_full_gradient_,
                        logger ? logger : this->logger_.get(), &myCpuTimeSec);
            have_cached_cost_ = true;
            have_cached_gradient_ = true;
        }
        // recycle old result
        std::copy(cached_gradient_.begin(), cached_gradient_.end(), gradient.begin());
        fval = cached_cost_;
    } else {
        if (!have_cached_cost_ || !std::equal(parameters.begin(), parameters.end(),
                                           cached_parameters_.begin())) {
            // Have to compute anew
            cached_status_ = hierarchical_wrapper_->evaluate(
                        parameters, cached_cost_, gsl::span<double>(),
                        cached_full_parameters_, cached_full_gradient_,
                        logger ? logger : this->logger_.get(), &myCpuTimeSec);
            have_cached_cost_ = true;
            have_cached_gradient_ = false;
        }
        fval = cached_cost_;
    }

    // update cached parameters
    cached_parameters_.resize(num_parameters_);
    std::copy(parameters.begin(), parameters.end(), cached_parameters_.begin());

    cpu_time_iteration_sec_ += myCpuTimeSec;
    cpu_time_total_sec_ += myCpuTimeSec;
    if(cpuTime)
        *cpuTime = myCpuTimeSec;

    if(afterCostFunctionCall(
                parameters, cached_cost_,
                gradient.data() ? cached_full_gradient_ : gsl::span<double>()
                ) != 0)
        return functionEvaluationFailure;

    return cached_status_;
}

void HierarchicalOptimizationReporter::finished(
        double optimalCost, gsl::span<const double> parameters, int exitStatus) const
{
    double timeElapsed = wall_timer_.getTotal();

    if(cached_cost_ > optimalCost) {
        // the optimal value is not from the cached parameters and we did not get
        // the optimal full parameter vector. since we don't know them, rather set to nan
        cached_full_parameters_.assign(cached_full_parameters_.size(), NAN);
        std::copy(parameters.begin(), parameters.end(), cached_parameters_.data());
        if(logger_) logger_->logmessage(LOGLVL_INFO, "cachedCost != optimalCost");
        cached_cost_ = NAN;
    }

    if(logger_)
        logger_->logmessage(LOGLVL_INFO, "Optimizer status %d, final llh: %e, time: wall: %f cpu: %f.",
                           exitStatus, cached_cost_, timeElapsed, cpu_time_total_sec_);

    if(result_writer_)
        result_writer_->saveOptimizerResults(cached_cost_, cached_full_parameters_,
                                           timeElapsed, cpu_time_total_sec_, exitStatus);
}

const std::vector<double> &HierarchicalOptimizationReporter::getFinalParameters() const
{
    return cached_full_parameters_;
}

bool HierarchicalOptimizationReporter::iterationFinished(
        gsl::span<const double> parameters,
        double objectiveFunctionValue,
        gsl::span<const double>  /*objectiveFunctionGradient*/) const
{
    double wallTimeIter = wall_timer_.getRound();
    double wallTimeOptim = wall_timer_.getTotal();

    if(logger_)
        logger_->logmessage(LOGLVL_INFO,
                           "iter: %d cost: %g "
                           "time_iter: wall: %gs cpu: %gs "
                           "time_optim: wall: %gs cpu: %gs",
                           num_iterations_, objectiveFunctionValue,
                           wallTimeIter, cpu_time_iteration_sec_,
                           wallTimeOptim, cpu_time_total_sec_);

    if(result_writer_) {
        /* check if the optimizer-reported cost matches the last function evaluation.
         * if so, we can log our cached parameter and gradient, otherwise we need to rely on what
         * the optimizer provided us. if no parameters are provided, we will still save the cached
         * one, even if the cost does not match, since this is the best parameter guess we have.
         */
        if(almostEqual(objectiveFunctionValue, cached_cost_)
                && (parameters.empty()
                    || std::equal(parameters.begin(), parameters.end(),
                                  cached_parameters_.begin()))) {
            result_writer_->logOptimizerIteration(
                        num_iterations_,
                        cached_full_parameters_,
                        objectiveFunctionValue,
                        // This might be misleading, the gradient could have been
                        // evaluated at other parameters if there was a line search inbetween
                        cached_full_gradient_,
                        wallTimeIter,
                        cpu_time_iteration_sec_);
        } else {
            // We don't have the full parameter vector, only the outer parameters
            // so we can't append them due to different dimension
            // TODO: save both, outer + combined? can easily save outer + inner separetly
            std::vector<double> nanParameters(cached_full_parameters_.size(), NAN);

            result_writer_->logOptimizerIteration(num_iterations_,
                                                nanParameters,
                                                objectiveFunctionValue,
                                                nanParameters,
                                                wallTimeIter,
                                                cpu_time_iteration_sec_);
        }
    }
    ++num_iterations_;

    logger_->setPrefix(default_logger_prefix_ + "i" + std::to_string(num_iterations_));
    cpu_time_iteration_sec_ = 0.0;

    return false;

}

bool HierarchicalOptimizationReporter::afterCostFunctionCall(
        gsl::span<const double>  /*parameters*/,
        double objectiveFunctionValue,
        gsl::span<const double> objectiveFunctionGradient) const
{
    double wallTime = wall_timer_.getTotal();
    //(double)(timeCostEvaluationEnd - timeCostEvaluationBegin) / CLOCKS_PER_SEC;

    if(!std::isfinite(objectiveFunctionValue))
        printObjectiveFunctionFailureMessage();

    if(result_writer_) {
        result_writer_->logObjectiveFunctionEvaluation(
                    cached_full_parameters_, cached_cost_,
                    objectiveFunctionGradient, num_iterations_,
                    num_function_calls_, wallTime);
    }
    return false;
}

void checkGradientForAnalyticalParameters(
        const std::vector<double> &gradient,
        const std::vector<int> &analyticalIndices, double threshold)
{
    for(auto const idx: analyticalIndices) {
        auto curGradient = gradient[idx];
        //std::cout<<"    : "<<idx<<"\t"<<curGradient<<std::endl;
        if(std::fabs(curGradient) > threshold)
            logmessage(LOGLVL_WARNING,
                       "Gradient w.r.t. analytically computed parameter "
                       "%d is %f, exceeding threshold %g", idx, curGradient,
                       threshold);
    }
}



} // namespace parpe
