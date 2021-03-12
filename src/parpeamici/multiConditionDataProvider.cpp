#include <parpeamici/multiConditionDataProvider.h>

#include <parpeamici/amiciMisc.h>
#include <parpecommon/logging.h>
#include <parpecommon/misc.h>
#include <parpecommon/parpeException.h>

#include <amici/amici.h>
#include <amici/hdf5.h>

#include <cassert>
#include <cmath>
#include <cstring>
#include <numeric>

namespace parpe {

MultiConditionDataProviderHDF5::MultiConditionDataProviderHDF5(
  std::unique_ptr<amici::Model> model,
  std::string const& hdf5Filename)
  : MultiConditionDataProviderHDF5(std::move(model), hdf5Filename, "")
{}

MultiConditionDataProviderHDF5::MultiConditionDataProviderHDF5(
  std::unique_ptr<amici::Model> model,
  std::string const& hdf5Filename,
  std::string const& rootPath)
  : model_(std::move(model))
  , root_path_(rootPath)
{

    [[maybe_unused]] auto lock = hdf5MutexGetLock();
    file_ = hdf5OpenForReading(hdf5Filename);

    optimization_options_ = parpe::OptimizationOptions::fromHDF5(file_);

    hdf5_measurement_path_ = rootPath + "/measurements/y";
    hdf5_measurement_sigma_path_ = rootPath + "/measurements/ysigma";
    hdf5_condition_path_ = rootPath + "/fixedParameters/k";
    hdf5_reference_condition_path_ =
      rootPath + "/fixedParameters/simulationConditions";
    hdf5_amici_options_path_ = rootPath + "/amiciOptions";
    hdf5_parameter_path_ = rootPath + "/parameters";
    hdf5_parameter_min_path_ = hdf5_parameter_path_ + "/lowerBound";
    hdf5_parameter_max_path_ = hdf5_parameter_path_ + "/upperBound";
    hdf5_parameter_scale_simulation_path_ =
      hdf5_parameter_path_ + "/pscaleSimulation";
    hdf5_parameter_scale_optimization_path_ =
      hdf5_parameter_path_ + "/pscaleOptimization";
    hdf5_simulation_to_optimization_parameter_mapping_path_ =
      rootPath + "/parameters/optimizationSimulationMapping";
    hdf5_parameter_overrides_path = rootPath + "/parameters/parameterOverrides";
    hdf5_parameter_ids_path_ = rootPath + "/parameters/parameterNames";
    hdf5_reinitialization_idxs_path_ =
        rootPath + "/fixedParameters/reinitializationIndices";
    checkDataIntegrity();

    amici::hdf5::readModelDataFromHDF5(
      file_, *this->model_, hdf5_amici_options_path_);
}

int
MultiConditionDataProviderHDF5::getNumberOfSimulationConditions() const
{
    // TODO: add additional layer for selection of condition indices (for
    // testing and later for minibatch)
    // -> won't need different file for testing/validation splits
    // TODO: cache

    [[maybe_unused]] auto lock = hdf5MutexGetLock();

    int d1, d2;
    hdf5GetDatasetDimensions(
      file_.getId(), hdf5_reference_condition_path_.c_str(), 2, &d1, &d2);

    return d1;
}

std::vector<int>
MultiConditionDataProviderHDF5::getSimulationToOptimizationParameterMapping(
  int conditionIdx) const
{
    std::string path = hdf5_simulation_to_optimization_parameter_mapping_path_;

    if (file_.nameExists(path)) {
        return hdf5Read2DIntegerHyperslab(
          file_, path, model_->np(), 1, 0, conditionIdx);
    }

    // return trivial default mapping
    std::vector<int> defaultMap(model_->np());
    std::iota(defaultMap.begin(), defaultMap.end(), 0);

    return defaultMap;
}

void
MultiConditionDataProviderHDF5::mapSimulationToOptimizationGradientAddMultiply(
  int conditionIdx,
  gsl::span<double const> simulation,
  gsl::span<double> optimization,
  gsl::span<double const> parameters,
  double coefficient) const
{
    auto mapping = getSimulationToOptimizationParameterMapping(conditionIdx);

    // Need to consider varying scaling
    auto scaleOpt = getParameterScaleOpt();
    auto scaleSim = getParameterScaleSim(conditionIdx);

    for (int i = 0; i < model_->np(); ++i) {
        // some model parameter are not mapped if there is no respective data
        if (mapping[i] >= 0) {
            double newGrad = applyChainRule(
              simulation[i], parameters[i], scaleSim[i], scaleOpt[mapping[i]]);
            optimization[mapping[i]] += coefficient * newGrad;
        }
    }
}

void
MultiConditionDataProviderHDF5::mapAndSetOptimizationToSimulationVariables(
  int conditionIdx,
  gsl::span<const double> optimization,
  gsl::span<double> simulation,
  gsl::span<amici::ParameterScaling> optimizationScale,
  gsl::span<amici::ParameterScaling> simulationScale) const
{
    auto mapping = getSimulationToOptimizationParameterMapping(conditionIdx);

    std::vector<double> overrides;
    if (file_.nameExists(hdf5_parameter_overrides_path)) {
        overrides.resize(model_->np());
        hdf5Read2DDoubleHyperslab(file_,
                                  hdf5_parameter_overrides_path,
                                  model_->np(),
                                  1,
                                  0,
                                  conditionIdx,
                                  overrides);
    }

    for (int i = 0; i < model_->np(); ++i) {
        if (mapping[i] >= 0) {
            // map from optimization parameters
            simulation[i] = getScaledParameter(
              getUnscaledParameter(optimization[mapping[i]],
                                   optimizationScale[mapping[i]]),
              simulationScale[i]);
        } else if (!overrides.empty()) {
            // TODO do we need to rescale here? or done in PEtab?
            simulation[i] = overrides[i];
        } else {
            simulation[i] = NAN;
        }
    }
}

std::vector<amici::ParameterScaling>
MultiConditionDataProviderHDF5::getParameterScaleOpt() const
{
    [[maybe_unused]] auto lock = hdf5MutexGetLock();
    auto resInt = amici::hdf5::getIntDataset1D(
      file_, hdf5_parameter_scale_optimization_path_);
    std::vector<amici::ParameterScaling> res(resInt.size());
    for (unsigned int i = 0; i < resInt.size(); ++i)
        res[i] = static_cast<amici::ParameterScaling>(resInt[i]);
    return res;
}

amici::ParameterScaling
MultiConditionDataProviderHDF5::getParameterScaleOpt(int parameterIdx) const
{
    auto res =
      hdf5Read1DIntegerHyperslab(
        file_, hdf5_parameter_scale_optimization_path_, 1, parameterIdx)
        .at(0);
    return static_cast<amici::ParameterScaling>(res);
}

std::vector<amici::ParameterScaling>
MultiConditionDataProviderHDF5::getParameterScaleSim(int simulationIdx) const
{
    auto resInt =
      hdf5Read2DIntegerHyperslab(file_,
                                 hdf5_parameter_scale_simulation_path_,
                                 1,
                                 model_->np(),
                                 simulationIdx,
                                 0);
    std::vector<amici::ParameterScaling> res(resInt.size());
    for (unsigned int i = 0; i < resInt.size(); ++i)
        res[i] = static_cast<amici::ParameterScaling>(resInt[i]);
    return res;
}

amici::ParameterScaling
MultiConditionDataProviderHDF5::getParameterScaleSim(
  int simulationIdx,
  int modelParameterIdx) const
{
    auto res = hdf5Read2DIntegerHyperslab(file_,
                                          hdf5_parameter_scale_simulation_path_,
                                          1,
                                          1,
                                          simulationIdx,
                                          modelParameterIdx)
                 .at(0);
    return static_cast<amici::ParameterScaling>(res);
}

void
MultiConditionDataProviderHDF5::updateFixedSimulationParameters(
  int simulationIdx,
  amici::ExpData& edata) const
{
    edata.fixedParameters.resize(model_->nk());

    // TODO cache
    int conditionIdxPreeq, conditionIdxSim;
    getSimAndPreeqConditions(simulationIdx,
                             conditionIdxPreeq,
                             conditionIdxSim);

    if (conditionIdxPreeq >= 0) {
        // -1 means no preequilibration
        edata.fixedParametersPreequilibration.resize(model_->nk());
        readFixedSimulationParameters(conditionIdxPreeq,
                                      edata.fixedParametersPreequilibration);
        edata.reinitialization_state_idxs_sim =
            getReinitializationIndices(simulationIdx);

    } else {
        edata.fixedParametersPreequilibration.resize(0);
        edata.reinitialization_state_idxs_sim .clear();
    }
    readFixedSimulationParameters(conditionIdxSim, edata.fixedParameters);
}

void MultiConditionDataProviderHDF5::setModel(std::unique_ptr<amici::Model> model)
{
    model_ = std::move(model);
}

std::vector<std::string> MultiConditionDataProviderHDF5::getProblemParameterIds() const
{
    return hdf5Read1dStringDataset(file_, hdf5_parameter_ids_path_);

}

void
MultiConditionDataProviderHDF5::readFixedSimulationParameters(
        int conditionIdx,
        gsl::span<double> buffer) const
{
    if (!model_->nk())
        return;

    [[maybe_unused]] auto lock = hdf5MutexGetLock();

    H5_SAVE_ERROR_HANDLER;

    hdf5Read2DDoubleHyperslab(file_.getId(),
                              hdf5_condition_path_.c_str(),
                              model_->nk(),
                              1,
                              0,
                              conditionIdx,
                              buffer);

    if (H5Eget_num(H5E_DEFAULT)) {
        logmessage(LOGLVL_CRITICAL,
                   "Problem in readFixedParameters (row %d, nk %d)\n",
                   conditionIdx,
                   model_->nk());
        printBacktrace(20);
        H5Ewalk2(
          H5E_DEFAULT, H5E_WALK_DOWNWARD, hdf5ErrorStackWalker_cb, nullptr);
        abort();
    }

    H5_RESTORE_ERROR_HANDLER;

    if (H5Eget_num(H5E_DEFAULT))
        throw ParPEException("MultiConditionDataProviderHDF5::"
                             "updateFixedSimulationParameters "
                             "unable to read data");
}

std::unique_ptr<amici::ExpData> MultiConditionDataProviderHDF5::getExperimentalDataForCondition(
        int simulationIdx) const {
    auto edata = std::make_unique<amici::ExpData>(*model_);

    [[maybe_unused]] auto lock = hdf5MutexGetLock();
    edata->setTimepoints(
                amici::hdf5::getDoubleDataset1D(
                    file_, root_path_ + "/measurements/t/"
                    + std::to_string(simulationIdx)));
    edata->setObservedData(getMeasurementForSimulationIndex(simulationIdx));
    edata->setObservedDataStdDev(getSigmaForSimulationIndex(simulationIdx));
    updateFixedSimulationParameters(simulationIdx, *edata);

    return edata;
}

std::vector<std::vector<double>>
MultiConditionDataProviderHDF5::getAllMeasurements() const
{
    std::vector<std::vector<double>> result(getNumberOfSimulationConditions());
    for (int conditionIdx = 0; (unsigned)conditionIdx < result.size();
         ++conditionIdx) {
        result[conditionIdx] = getMeasurementForSimulationIndex(conditionIdx);
    }
    return result;
}

std::vector<std::vector<double>>
MultiConditionDataProviderHDF5::getAllSigmas() const
{
    std::vector<std::vector<double>> result(getNumberOfSimulationConditions());
    for (int conditionIdx = 0; (unsigned)conditionIdx < result.size();
         ++conditionIdx) {
        result[conditionIdx] = getSigmaForSimulationIndex(conditionIdx);
    }
    return result;
}

std::vector<double>
MultiConditionDataProviderHDF5::getSigmaForSimulationIndex(
  int simulationIdx) const
{
    hsize_t dim1, dim2;
    [[maybe_unused]] auto lock = hdf5MutexGetLock();
    return amici::hdf5::getDoubleDataset2D(file_,
                                           hdf5_measurement_sigma_path_ + "/" +
                                             std::to_string(simulationIdx),
                                           dim1,
                                           dim2);
}

std::vector<double>
MultiConditionDataProviderHDF5::getMeasurementForSimulationIndex(
  int simulationIdx) const
{
    hsize_t dim1, dim2;
    [[maybe_unused]] auto lock = hdf5MutexGetLock();
    return amici::hdf5::getDoubleDataset2D(file_,
                                           hdf5_measurement_path_ + "/" +
                                             std::to_string(simulationIdx),
                                           dim1,
                                           dim2);
}

void
MultiConditionDataProviderHDF5::getOptimizationParametersLowerBounds(
  gsl::span<double> buffer) const
{
    [[maybe_unused]] auto lock = hdf5MutexGetLock();

    auto dataset = file_.openDataSet(hdf5_parameter_min_path_);

    auto dataspace = dataset.getSpace();
    // hdf5ParameterMinPath dimensions don't match
    Expects(dataspace.getSimpleExtentNdims() == 1);
    hsize_t dim = 0;
    dataspace.getSimpleExtentDims(&dim);
    Expects(dim == (unsigned)getNumOptimizationParameters());
    Expects(dim == buffer.size());
    dataset.read(buffer.data(), H5::PredType::NATIVE_DOUBLE);
}

void
MultiConditionDataProviderHDF5::getOptimizationParametersUpperBounds(
  gsl::span<double> buffer) const
{
    [[maybe_unused]] auto lock = hdf5MutexGetLock();

    auto dataset = file_.openDataSet(hdf5_parameter_max_path_);

    auto dataspace = dataset.getSpace();
    // hdf5ParameterMaxPath dimensions dont match
    Expects(dataspace.getSimpleExtentNdims() == 1);
    hsize_t dim = 0;
    dataspace.getSimpleExtentDims(&dim);
    Expects(dim == (unsigned)getNumOptimizationParameters());
    Expects(dim == buffer.size());
    dataset.read(buffer.data(), H5::PredType::NATIVE_DOUBLE);
}

int
MultiConditionDataProviderHDF5::getNumOptimizationParameters() const
{
    std::string path = root_path_ + "/parameters/parameterNames";
    int size = 0;
    hdf5GetDatasetDimensions(file_.getId(), path.c_str(), 1, &size);
    return size;
}

std::unique_ptr<amici::Model>
MultiConditionDataProviderHDF5::getModel() const
{
    return std::unique_ptr<amici::Model>(model_->clone());
}

std::unique_ptr<amici::Solver>
MultiConditionDataProviderHDF5::getSolver() const
{
    auto solver = model_->getSolver();
    [[maybe_unused]] auto lock = hdf5MutexGetLock();

    amici::hdf5::readSolverSettingsFromHDF5(
      file_, *solver, hdf5_amici_options_path_);
    return solver;
}

void
MultiConditionDataProviderHDF5::updateSimulationParametersAndScale(
  int simulationIdx,
  gsl::span<const double> optimizationParams,
  amici::Model& model) const
{
    // int conditionIdxPreeq, conditionIdxSim;
    // getSimAndPreeqConditions(simulationIdx, conditionIdxPreeq,
    // conditionIdxSim);

    auto scaleSim = getParameterScaleSim(simulationIdx);
    auto p = model.getParameters();
    auto scaleOpt = getParameterScaleOpt();

    model.setParameterScale(scaleSim);
    mapAndSetOptimizationToSimulationVariables(
      simulationIdx, optimizationParams, p, scaleOpt, scaleSim);
    model.setParameters(p);
}

void
MultiConditionDataProviderHDF5::copyInputData(H5::H5File const& target)
{

    H5Ocopy(file_.getId(),
            "/",
            target.getId(),
            "/inputData",
            H5P_DEFAULT,
            H5P_DEFAULT);
    H5Fflush(target.getId(), H5F_SCOPE_LOCAL);
}

void
MultiConditionDataProviderHDF5::getSimAndPreeqConditions(
  const int simulationIdx,
  int& preequilibrationConditionIdx,
  int& simulationConditionIdx) const
{
    auto tmp = hdf5Read2DIntegerHyperslab(
      file_, hdf5_reference_condition_path_, 1, 3, simulationIdx, 0);
    preequilibrationConditionIdx = tmp[0];
    simulationConditionIdx = tmp[1];
}

std::vector<int> MultiConditionDataProviderHDF5::getReinitializationIndices(
    const int simulationIdx) const {
    [[maybe_unused]] auto lock = hdf5MutexGetLock();
    auto dataset = file_.openDataSet(hdf5_reinitialization_idxs_path_);

    auto filespace = dataset.getSpace();
    Expects(filespace.getSimpleExtentNdims() == 1);
    hsize_t num_simulation_conditions;
    filespace.getSimpleExtentDims(&num_simulation_conditions);
    Expects(simulationIdx >= 0
            && (hsize_t) simulationIdx < num_simulation_conditions);

    // read only for one condition
    const hsize_t len = 1;
    const hsize_t offset = simulationIdx;
    filespace.selectHyperslab(H5S_SELECT_SET, &len, &offset);
    H5::DataSpace memspace(1, &len);

    auto memtype = H5::VarLenType(H5::PredType::NATIVE_INT);

    hvl_t buffer;
    dataset.read(&buffer, memtype, memspace, filespace);

    Expects(buffer.p);
    auto int_ptr  = static_cast<int *>(buffer.p);

    return std::vector<int>(&int_ptr[0], &int_ptr[buffer.len]);
}

H5::H5File MultiConditionDataProviderHDF5::getHdf5File() const
{
    [[maybe_unused]] auto lock = hdf5MutexGetLock();
    H5::H5File result(file_);
    return result;
}

// void MultiConditionDataProvider::printInfo() const {
//    int maxwidth = 25;
//    int numMultiStartRuns = getNumMultiStartRuns();
//    logmessage(LOGLVL_INFO, "%*s: %d", maxwidth, "Num multistart optims",
//    numMultiStartRuns);

//    for(int i = 0; i < numMultiStartRuns; ++i) {
//        int numStarts = getNumLocalOptimizationsForMultiStartRun(i);
//        logmessage(LOGLVL_INFO, "%*s: %d", maxwidth, "Num starts",
//        numStarts);
//        logmessage(LOGLVL_INFO, "%*s: %d", maxwidth, "Genotypes",
//        getNumGenotypes());
//    }

//    logmessage(LOGLVL_INFO, "%*s: %d", maxwidth, "Max iterations",
//    getMaxIter());
//}

void
MultiConditionDataProviderHDF5::checkDataIntegrity() const
{
    // check matching IDs
    std::string modelParameterIdsPath = root_path_ + "/model/parameterIds";
    auto dataParameterIds =
      hdf5Read1dStringDataset(file_, modelParameterIdsPath);
    auto modelParameterIds = this->model_->getParameterIds();
    RELEASE_ASSERT(dataParameterIds == modelParameterIds,
                   "Parameter IDs do not match.");

    std::string speciesIdsPath = root_path_ + "/model/stateIds";
    auto dataSpeciesIds = hdf5Read1dStringDataset(file_, speciesIdsPath);
    RELEASE_ASSERT(dataSpeciesIds == model_->getStateIds(),
                   "State IDs do not match.");

    std::string fixedParIdsPath = root_path_ + "/model/fixedParameterIds";
    auto fixedParIds = hdf5Read1dStringDataset(file_, fixedParIdsPath);
    RELEASE_ASSERT(fixedParIds == model_->getFixedParameterIds(),
                   "Fixed parameter IDs do not match.");

    std::string observableIdsPath = root_path_ + "/model/observableIds";
    auto observableIds = hdf5Read1dStringDataset(file_, observableIdsPath);
    RELEASE_ASSERT(observableIds == model_->getObservableIds(),
                   "Observable IDs do not match.");

    // int numConditions = getNumberOfSimulationConditions();

    int d1, d2; //, d3;

    [[maybe_unused]] auto lock = hdf5MutexGetLock();

    Ensures(
      H5Lexists(file_.getId(), hdf5_measurement_path_.c_str(), H5P_DEFAULT));
    Ensures(H5Lexists(
      file_.getId(), hdf5_measurement_sigma_path_.c_str(), H5P_DEFAULT));

    //    parpe::hdf5GetDatasetDimensions(file.getId(),
    //    hdf5MeasurementPath.c_str(),
    //                                    3, &d1, &d2, &d3);
    //    RELEASE_ASSERT(d1 >= numConditions, "");
    //    RELEASE_ASSERT(d2 == model->nytrue, "");
    //    RELEASE_ASSERT(d3 >= model->nt(), "");

    //    parpe::hdf5GetDatasetDimensions(file.getId(),
    //    hdf5MeasurementSigmaPath.c_str(),
    //                                    3, &d1, &d2, &d3);
    //    RELEASE_ASSERT(d1 >= numConditions, "");
    //    RELEASE_ASSERT(d2 == model->nytrue, "");
    //    RELEASE_ASSERT(d3 >= model->nt(), "");

    if (model_->nk()) {
        parpe::hdf5GetDatasetDimensions(
          file_.getId(), hdf5_condition_path_.c_str(), 2, &d1, &d2);
        Expects(d1 == model_->nk());
    }
}

MultiConditionDataProviderDefault::MultiConditionDataProviderDefault(
  std::unique_ptr<amici::Model> model,
  std::unique_ptr<amici::Solver> solver)
  : model_(std::move(model))
  , solver_(std::move(solver))
{}

int
MultiConditionDataProviderDefault::getNumberOfSimulationConditions() const
{
    return edata_.size();
}

std::vector<int>
MultiConditionDataProviderDefault::getSimulationToOptimizationParameterMapping(
  int /*conditionIdx*/) const
{
    std::vector<int> mapping(model_->np());
    std::iota(mapping.begin(), mapping.end(), 0);
    return mapping;
}

void
MultiConditionDataProviderDefault ::
  mapSimulationToOptimizationGradientAddMultiply(
    int conditionIdx,
    gsl::span<double const> simulation,
    gsl::span<double> optimization,
    gsl::span<const double> /*parameters*/,
    double coefficient) const
{
    // TODO redundant
    auto mapping = getSimulationToOptimizationParameterMapping(conditionIdx);

    for (int i = 0; i < model_->np(); ++i) {
        optimization[mapping[i]] = coefficient * simulation[i];
    }
}

void
MultiConditionDataProviderDefault::mapAndSetOptimizationToSimulationVariables(
  int conditionIdx,
  gsl::span<const double> optimization,
  gsl::span<double> simulation,
  gsl::span<amici::ParameterScaling> /*optimizationScale*/,
  gsl::span<amici::ParameterScaling> /*simulationScale*/) const
{
    // TODO redundant
    auto mapping = getSimulationToOptimizationParameterMapping(conditionIdx);

    for (int i = 0; i < model_->np(); ++i) {
        simulation[i] = optimization[mapping[i]];
    }
}

std::vector<amici::ParameterScaling>
MultiConditionDataProviderDefault::getParameterScaleOpt() const
{
    return model_->getParameterScale();
}

amici::ParameterScaling
MultiConditionDataProviderDefault::getParameterScaleOpt(
  int optimizationParameterIndex) const
{
    return getParameterScaleSim(0, optimizationParameterIndex);
}

amici::ParameterScaling
MultiConditionDataProviderDefault::getParameterScaleSim(
  int /*simulationIdx*/,
  int optimizationParameterIndex) const
{
    // TODO assumes no extra optimization parameters
    return model_->getParameterScale()[optimizationParameterIndex];
}

std::vector<amici::ParameterScaling>
MultiConditionDataProviderDefault::getParameterScaleSim(
  int /*simulationIdx*/) const
{
    return model_->getParameterScale();
}

void
MultiConditionDataProviderDefault::updateSimulationParametersAndScale(
  int /*conditionIndex*/,
  gsl::span<const double> optimizationParams,
  amici::Model& model) const
{
    logmessage(LOGLVL_WARNING,
               "MultiConditionDataProviderDefault::updateSimulationParameters: "
               "No proper mapping implemented. Ensure this is correct.");
    model.setParameters(std::vector<double>(optimizationParams.begin(),
                                            optimizationParams.end()));
}

std::unique_ptr<amici::ExpData>
MultiConditionDataProviderDefault::getExperimentalDataForCondition(
  int conditionIdx) const
{
    return std::make_unique<amici::ExpData>(edata_[conditionIdx]);
}

std::vector<std::vector<double>>
MultiConditionDataProviderDefault::getAllMeasurements() const
{
    std::vector<std::vector<double>> measurements;
    for (const auto& e : edata_) {
        measurements.push_back(e.getObservedData());
    }
    return measurements;
}

std::vector<std::vector<double>>
MultiConditionDataProviderDefault::getAllSigmas() const
{
    std::vector<std::vector<double>> sigmas;
    for (const auto& e : edata_) {
        sigmas.push_back(e.getObservedDataStdDev());
    }
    return sigmas;
}

int
MultiConditionDataProviderDefault::getNumOptimizationParameters() const
{
    return model_->np();
}

std::unique_ptr<amici::Model>
MultiConditionDataProviderDefault::getModel() const
{
    return std::unique_ptr<amici::Model>(model_->clone());
}

std::unique_ptr<amici::Solver>
MultiConditionDataProviderDefault::getSolver() const
{
    return std::unique_ptr<amici::Solver>(solver_->clone());
}

std::vector<std::string> MultiConditionDataProviderDefault::getProblemParameterIds() const
{
    // not implemented
    std::terminate();
}

double
applyChainRule(double gradient,
               double parameter,
               amici::ParameterScaling oldScale,
               amici::ParameterScaling newScale)
{
    if (oldScale == newScale)
        return gradient;

    // unapply old
    switch (oldScale) {
        case amici::ParameterScaling::log10:
            gradient /= getUnscaledParameter(parameter, oldScale) * log(10);
            break;
        case amici::ParameterScaling::ln:
            gradient /= getUnscaledParameter(parameter, oldScale);
            break;
        case amici::ParameterScaling::none:
            break;
    }

    // apply
    switch (newScale) {
        case amici::ParameterScaling::log10:
            gradient *= getUnscaledParameter(parameter, oldScale) * log(10);
            break;
        case amici::ParameterScaling::ln:
            gradient *= getUnscaledParameter(parameter, oldScale);
            break;
        case amici::ParameterScaling::none:
            break;
    }

    return gradient;
}

} // namespace parpe
