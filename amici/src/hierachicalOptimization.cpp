#include "hierachicalOptimization.h"

#ifdef __INTEL_COMPILER
// constexpr did not work on icc (ICC) 16.0.4 20160811
#define constexpr
#endif

namespace parpe {


HierachicalOptimizationWrapper::HierachicalOptimizationWrapper(std::unique_ptr<AmiciSummedGradientFunction<int> > fun, std::unique_ptr<parpe::ScalingFactorHdf5Reader> reader, int numConditions, int numObservables, int numTimepoints)
    : fun(fun.release()),
      reader(std::move(reader)),
      numConditions(numConditions),
      numObservables(numObservables),
      numTimepoints(numTimepoints) {

    proportionalityFactorIndices = this->reader->getProportionalityFactorIndices();
    std::sort(this->proportionalityFactorIndices.begin(),
              this->proportionalityFactorIndices.end());
}

FunctionEvaluationStatus HierachicalOptimizationWrapper::evaluate(const double * const parameters, double &fval, double *gradient) const {

    if(numScalingFactors() == 0) {
        // nothing to do, just pass through
        std::vector<int> dataIndices(numConditions);
        std::iota(dataIndices.begin(), dataIndices.end(), 0);

        return fun->evaluate(parameters, dataIndices, fval, gradient);
    }

    // evaluate with scaling parameters set to 1
    auto modelOutput = getModelOutputs(parameters);

    // compute correct scaling factors analytically
    auto scalings = computeAnalyticalScalings(modelOutput);

    // evaluate with analytical scaling parameters
    auto status = evaluateWithScalings(parameters, scalings, modelOutput, fval, gradient);

    return status;
}

std::vector<std::vector<double> > HierachicalOptimizationWrapper::getModelOutputs(const double * const reducedParameters) const {
    // run simulations (no gradient!) with scaling parameters == 1, collect outputs

    std::vector<double> scalingDummy(numScalingFactors(), 1);
    auto fullParameters = getFullParameters(reducedParameters, scalingDummy);

    std::vector<std::vector<double> > modelOutput(numConditions);
    fun->getModelOutputs(fullParameters.data(), modelOutput);

    return modelOutput;
}

std::vector<double> HierachicalOptimizationWrapper::computeAnalyticalScalings(std::vector<std::vector<double> > &modelOutputs) const {
    // NOTE: does not handle replicates, assumes normal distribution, does not compute sigmas

    int numProportionalityFactors = proportionalityFactorIndices.size();
    std::vector<double> proportionalityFactors(numProportionalityFactors);

    auto measurements = fun->getAllMeasurements();

    for(int i = 0; i < numProportionalityFactors; ++i) {
        proportionalityFactors[i] = optimalScaling(i, modelOutputs, measurements);
        applyOptimalScaling(i, proportionalityFactors[i], modelOutputs);
    }

    return proportionalityFactors;
}

void HierachicalOptimizationWrapper::applyOptimalScaling(int scalingIdx, double scaling, std::vector<std::vector<double> > &modelOutputs) const {
    auto dependentConditions = reader->getConditionsForScalingParameter(scalingIdx);
    for (auto const conditionIdx: dependentConditions) {
        auto dependentObservables = reader->getObservablesForScalingParameter(scalingIdx, conditionIdx);
        for(auto const observableIdx: dependentObservables) {
            for(int timeIdx = 0; timeIdx < numTimepoints; ++timeIdx) {
                modelOutputs[conditionIdx][observableIdx + timeIdx * numObservables] *= scaling;
            }
        }
    }
}

double HierachicalOptimizationWrapper::optimalScaling(int scalingIdx, const std::vector<std::vector<double> > &modelOutputsUnscaled, const std::vector<std::vector<double> > &measurements) const {
    auto dependentConditions = reader->getConditionsForScalingParameter(scalingIdx);

    double enumerator = 0.0;
    double denominator = 0.0;

    for (auto const conditionIdx: dependentConditions) {
        auto dependentObservables = reader->getObservablesForScalingParameter(scalingIdx, conditionIdx);
        for(auto const observableIdx: dependentObservables) {
            for(int timeIdx = 0; timeIdx < numTimepoints; ++timeIdx) {
                double mes = measurements[conditionIdx][observableIdx + timeIdx * numObservables];
                if(!std::isnan(mes)) {
                    double sim = modelOutputsUnscaled[conditionIdx][observableIdx + timeIdx * numObservables];
                    enumerator += sim * mes;
                    denominator += sim * sim;
                }
            }
        }
    }

    double proportionalityFactor = enumerator / denominator;

    return proportionalityFactor;
}

FunctionEvaluationStatus HierachicalOptimizationWrapper::evaluateWithScalings(const double * const reducedParameters, const std::vector<double> &scalings, const std::vector<std::vector<double> > &modelOutputsScaled, double &fval, double *gradient) const {

    if(gradient) {
        // simulate with updated theta for sensitivities
        std::fill(gradient, gradient + numParameters(), 0.0);
        auto fullParameters = getFullParameters(reducedParameters, scalings);
        // TODO: only those necessary
        std::vector<int> dataIndices(numConditions);
        std::iota(dataIndices.begin(), dataIndices.end(), 0);
        fun->evaluate(fullParameters.data(), dataIndices, fval, gradient);

    } else {
        // just compute
        fval = computeLikelihood(modelOutputsScaled);
    }

    return functionEvaluationSuccess;
}

double HierachicalOptimizationWrapper::computeLikelihood(const std::vector<std::vector<double> > &modelOutputsScaled) const {
    double llh = 0.0;
    constexpr double pi = atan(1)*4;
    auto measurements = fun->getAllMeasurements();
    double sigmaSquared = 1.0; // NOTE: no user-provided sigma supported at the moment

    for (int conditionIdx = 0; conditionIdx < numConditions; ++conditionIdx) {
        for(int observableIdx = 0; observableIdx < numObservables; ++observableIdx) {
            for(int timeIdx = 0; timeIdx < numTimepoints; ++timeIdx) {
                double mes = measurements[conditionIdx][observableIdx + timeIdx * numObservables];
                if(!std::isnan(mes)) {
                    double sim = modelOutputsScaled[conditionIdx][observableIdx + timeIdx * numObservables];
                    double diff = mes - sim;
                    diff *= diff;
                    llh += log(2.0 * pi * sigmaSquared) + diff / sigmaSquared;
                }
            }
        }
    }

    llh /= 2.0;

    return llh;
}

std::vector<double> HierachicalOptimizationWrapper::getFullParameters(const double * const reducedParameters, const std::vector<double> &scalingFactors) const {
    std::vector<double> fullParameters(fun->numParameters());
    int idxScaling = 0;
    int idxRegular = 0;

    for(int i = 0; i < (signed) fullParameters.size(); ++i) {
        if(proportionalityFactorIndices[idxScaling] == i)
            fullParameters[i] = scalingFactors[idxScaling++];
        else
            fullParameters[i] = reducedParameters[idxRegular++];
    }

    return fullParameters;
}

int HierachicalOptimizationWrapper::numParameters() const {
    return fun->numParameters() - numScalingFactors();
}

int HierachicalOptimizationWrapper::numScalingFactors() const {
    return proportionalityFactorIndices.size();
}

const std::vector<int> &HierachicalOptimizationWrapper::getProportionalityFactorIndices() const
{
    return proportionalityFactorIndices;
}

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

    return proportionalityFactorIndices;
}

void ScalingFactorHdf5Reader::readFromFile() {
    auto lock = hdf5MutexGetLock();
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

}

HierachicalOptimizationProblemWrapper::HierachicalOptimizationProblemWrapper(std::unique_ptr<OptimizationProblem> problemToWrap, MultiConditionDataProvider *dataProvider)
    : wrappedProblem(std::move(problemToWrap))
{
    auto wrappedFun = static_cast<SummedGradientFunctionGradientFunctionAdapter<int>*>(wrappedProblem->costFun.get());
    costFun.reset(new HierachicalOptimizationWrapper(
                      std::unique_ptr<AmiciSummedGradientFunction<int>>(
                          static_cast<AmiciSummedGradientFunction<int>*>(wrappedFun->gradFun.get())),
                      std::make_unique<ScalingFactorHdf5Reader>(dataProvider->file, "/"),
                      dataProvider->getNumberOfConditions(),
                      dataProvider->model->nytrue,
                      dataProvider->model->nt()));
}

HierachicalOptimizationProblemWrapper::~HierachicalOptimizationProblemWrapper()
{
    costFun.release(); // Avoid double delete. This will be destroyed when wrappedProblem goes out of scope!
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
    HierachicalOptimizationWrapper* hierarchical = static_cast<HierachicalOptimizationWrapper *>(wrappedProblem->costFun.get());
    auto proportionalityFactorIndices = hierarchical->getProportionalityFactorIndices();
    unsigned int nextScalingIdx = 0;
    unsigned int bufferIdx = 0;
    for(int i = 0; (unsigned)i < fullParams.size(); ++i) {
        if(proportionalityFactorIndices[nextScalingIdx] == i) {
            // skip
            ++nextScalingIdx;
        } else {
            // copy
            buffer[bufferIdx] = fullParams[i];
            ++bufferIdx;
        }
    }
    assert(nextScalingIdx == proportionalityFactorIndices.size());
    assert(bufferIdx == (unsigned)costFun->numParameters());
}


} // namespace parpe
