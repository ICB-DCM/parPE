#include <parpeamici/multiConditionDataProvider.h>

#include <parpecommon/logging.h>
#include <parpecommon/misc.h>
#include <parpecommon/parpeException.h>

#include <amici/amici.h>
#include <amici/hdf5.h>

#include <cassert>
#include <cstring>
#include <cmath>
#include <numeric>

namespace parpe {

MultiConditionDataProviderHDF5::MultiConditionDataProviderHDF5(
        std::unique_ptr<amici::Model> model,
        std::string const& hdf5Filename)
    : MultiConditionDataProviderHDF5(std::move(model), hdf5Filename, "") {}

MultiConditionDataProviderHDF5::MultiConditionDataProviderHDF5(
        std::unique_ptr<amici::Model> model,
        std::string const& hdf5Filename,
        std::string const& rootPath)
    : model(std::move(model)), rootPath(rootPath) {

    auto lock = hdf5MutexGetLock();
    file = hdf5OpenForReading(hdf5Filename);

    optimizationOptions = parpe::OptimizationOptions::fromHDF5(getHdf5FileId());

    hdf5MeasurementPath = rootPath + "/measurements/y";
    hdf5MeasurementSigmaPath = rootPath + "/measurements/ysigma";
    hdf5ConditionPath = rootPath + "/fixedParameters/k";
    hdf5ReferenceConditionPath = rootPath
            + "/fixedParameters/simulationConditions";
    hdf5AmiciOptionPath = rootPath + "/amiciOptions";
    hdf5ParameterPath = rootPath + "/parameters";
    hdf5ParameterMinPath = hdf5ParameterPath + "/lowerBound";
    hdf5ParameterMaxPath = hdf5ParameterPath + "/upperBound";
    hdf5ParameterScalingPath = hdf5ParameterPath + "/pscale";
    hdf5SimulationToOptimizationParameterMappingPath = rootPath
            + "/parameters/optimizationSimulationMapping";

    amici::hdf5::readModelDataFromHDF5(file, *this->model, hdf5AmiciOptionPath);
}

int MultiConditionDataProviderHDF5::getNumberOfSimulationConditions() const {
    // TODO: add additional layer for selection of condition indices (for testing
    // and later for minibatch)
    // -> won't need different file for testing/validation splits
    // TODO: cache

    auto lock = hdf5MutexGetLock();

    int d1, d2;
    hdf5GetDatasetDimensions(file.getId(), hdf5ReferenceConditionPath.c_str(),
                             2, &d1, &d2);

    return d1;
}



std::vector<int> MultiConditionDataProviderHDF5::getSimulationToOptimizationParameterMapping(int conditionIdx) const  {
    std::string path = hdf5SimulationToOptimizationParameterMappingPath;

    if(hdf5DatasetExists(file.getId(), path)) {
        return hdf5Read2DIntegerHyperslab(file, path, model->np(), 1, 0, conditionIdx);
    }

    // return trivial default mapping
    std::vector<int> defaultMap(model->np());
    std::iota(defaultMap.begin(), defaultMap.end(), 0);

    return defaultMap;
}

void MultiConditionDataProviderHDF5::mapSimulationToOptimizationVariablesAddMultiply(
        int conditionIdx, gsl::span<double const> simulation,
        gsl::span<double> optimization, double coefficient) const {
    auto mapping = getSimulationToOptimizationParameterMapping(conditionIdx);

    for(int i = 0; i < model->np(); ++i) {
        // some model parameter are not mapped if there is no respective data
        if(mapping[i] >= 0)
            optimization[mapping[i]] += coefficient * simulation[i];
        else
            RELEASE_ASSERT(std::isnan(simulation[i]) || simulation[i] == 0.0, "");
    }
}

void MultiConditionDataProviderHDF5::mapAndSetOptimizationToSimulationVariables(
        int conditionIdx, gsl::span<const double> optimization,
        gsl::span<double> simulation) const
{
    auto mapping = getSimulationToOptimizationParameterMapping(conditionIdx);

    for(int i = 0; i < model->np(); ++i) {
        if(mapping[i] >= 0)
            simulation[i] = optimization[mapping[i]];
        else
            simulation[i] = NAN;
    }
}

amici::ParameterScaling MultiConditionDataProviderHDF5::getParameterScale(
        int optimizationParameterIndex) const
{
    auto res = hdf5Read1DIntegerHyperslab(file, hdf5ParameterScalingPath,
                                          1, optimizationParameterIndex).at(0);
    return static_cast<amici::ParameterScaling>(res);
}

/**
 * @brief Update the contstants in AMICI ExpData object. Reads a slab for the
 * given simulation from fixed parameters matrix.
 *
 * @param simulationIdx Index of the experimental condition for which the
 * parameters should be taken.
 * @param edata The object to be updated.
 */
void MultiConditionDataProviderHDF5::updateFixedSimulationParameters(
        int simulationIdx, amici::ExpData &edata) const {
    edata.fixedParameters.resize(model->nk());

    // TODO cache
    auto tmp = hdf5Read2DIntegerHyperslab(file.getId(),
                                          hdf5ReferenceConditionPath,
                                          1, 2, simulationIdx, 0);
    int conditionIdxPreeq = tmp[0];
    int conditionIdxSim = tmp[1];
    if(conditionIdxPreeq >= 0) {
        // -1 means no preequilibration
        edata.fixedParametersPreequilibration.resize(model->nk());
        readFixedSimulationParameters(
                    conditionIdxPreeq,
                    edata.fixedParametersPreequilibration.data());
    }
    readFixedSimulationParameters(conditionIdxSim,
                                  edata.fixedParameters.data());
}

void MultiConditionDataProviderHDF5::readFixedSimulationParameters(
        int conditionIdx, double *buffer) const
{
    if(!model->nk())
        return;

    auto lock = hdf5MutexGetLock();

    H5_SAVE_ERROR_HANDLER;

    hdf5Read2DDoubleHyperslab(file.getId(), hdf5ConditionPath.c_str(),
                              model->nk(), 1,
                              0, conditionIdx, buffer);

    if (H5Eget_num(H5E_DEFAULT)) {
        logmessage(LOGLVL_CRITICAL,
                   "Problem in readFixedParameters (row %d, nk %d)\n",
                   conditionIdx, model->nk());
        printBacktrace(20);
        H5Ewalk2(H5E_DEFAULT, H5E_WALK_DOWNWARD, hdf5ErrorStackWalker_cb, nullptr);
        abort();
    }

    H5_RESTORE_ERROR_HANDLER;

    if(H5Eget_num(H5E_DEFAULT))
        throw ParPEException("MultiConditionDataProviderHDF5::updateFixedSimulationParameters unable to read data");
}

std::unique_ptr<amici::ExpData> MultiConditionDataProviderHDF5::getExperimentalDataForCondition(
        int simulationIdx) const {
    auto lock = hdf5MutexGetLock();

    auto edata = std::make_unique<amici::ExpData>(*model);
    RELEASE_ASSERT(edata, "Failed getting experimental data. Check data file.");
    edata->setTimepoints(
                amici::hdf5::getDoubleDataset1D(
                    file.getId(), rootPath + "/measurements/t/"
                    + std::to_string(simulationIdx)));
    edata->setObservedData(getMeasurementForSimulationIndex(simulationIdx));
    edata->setObservedDataStdDev(getSigmaForSimulationIndex(simulationIdx));
    updateFixedSimulationParameters(simulationIdx, *edata);

    return edata;
}

std::vector<std::vector<double> > MultiConditionDataProviderHDF5::getAllMeasurements() const {
    std::vector<std::vector<double>> result(getNumberOfSimulationConditions());
    for(int conditionIdx = 0; (unsigned) conditionIdx < result.size(); ++conditionIdx) {
        result[conditionIdx] = getMeasurementForSimulationIndex(conditionIdx);
    }
    return result;
}

std::vector<std::vector<double> > MultiConditionDataProviderHDF5::getAllSigmas() const
{
    // TODO: how to deal with sigma parameters vs table
    std::vector<std::vector<double>> result(getNumberOfSimulationConditions());
    for(int conditionIdx = 0; (unsigned) conditionIdx < result.size(); ++conditionIdx) {
        result[conditionIdx]= getSigmaForSimulationIndex(conditionIdx);
    }
    return result;
}

std::vector<double> MultiConditionDataProviderHDF5::getSigmaForSimulationIndex(int simulationIdx) const
{
    hsize_t dim1, dim2;
    return amici::hdf5::getDoubleDataset2D(
                file.getId(),
                hdf5MeasurementSigmaPath + "/" + std::to_string(simulationIdx),
                dim1, dim2);
}

std::vector<double> MultiConditionDataProviderHDF5::getMeasurementForSimulationIndex(int simulationIdx) const
{
    hsize_t dim1, dim2;
    return amici::hdf5::getDoubleDataset2D(
                file.getId(),
                hdf5MeasurementPath + "/" + std::to_string(simulationIdx),
                dim1, dim2);
}

void MultiConditionDataProviderHDF5::getOptimizationParametersLowerBounds(
        double *buffer) const {
    auto lock = hdf5MutexGetLock();

    auto dataset = file.openDataSet(hdf5ParameterMinPath);

    auto dataspace = dataset.getSpace();
    RELEASE_ASSERT(dataspace.getSimpleExtentNdims() == 1,
                   "hdf5ParameterMinPath dimensions dont match");
    hsize_t dim = 0;
    dataspace.getSimpleExtentDims(&dim);
    RELEASE_ASSERT(dim == (unsigned) getNumOptimizationParameters(),
                   "hdf5ParameterMinPath dimensions dont match");

    dataset.read(buffer, H5::PredType::NATIVE_DOUBLE);
}

void MultiConditionDataProviderHDF5::getOptimizationParametersUpperBounds(
        double *buffer) const {
    auto lock = hdf5MutexGetLock();

    auto dataset = file.openDataSet(hdf5ParameterMaxPath);

    auto dataspace = dataset.getSpace();
    RELEASE_ASSERT(dataspace.getSimpleExtentNdims() == 1,
                   "hdf5ParameterMaxPath dimensions dont match");
    hsize_t dim = 0;
    dataspace.getSimpleExtentDims(&dim);
    RELEASE_ASSERT(dim == (unsigned) getNumOptimizationParameters(),
                   "hdf5ParameterMaxPath dimensions dont match");

    dataset.read(buffer, H5::PredType::NATIVE_DOUBLE);
}

int MultiConditionDataProviderHDF5::getNumOptimizationParameters() const {
    std::string path = rootPath + "/parameters/parameterNames";
    int size = 0;
    hdf5GetDatasetDimensions(file.getId(), path.c_str(), 1, &size);
    return size;
}


std::unique_ptr<amici::Model> MultiConditionDataProviderHDF5::getModel() const {
    return std::unique_ptr<amici::Model>(model->clone());
}

std::unique_ptr<amici::Solver> MultiConditionDataProviderHDF5::getSolver() const
{
    auto solver = model->getSolver();
    auto lock = hdf5MutexGetLock();

    amici::hdf5::readSolverSettingsFromHDF5(file.getId(), *solver, hdf5AmiciOptionPath);
    return solver;
}


void MultiConditionDataProviderHDF5::updateSimulationParameters(
        int conditionIndex, gsl::span<const double> optimizationParams,
        amici::Model &model) const
{
    auto p = model.getParameters();
    mapAndSetOptimizationToSimulationVariables(conditionIndex, optimizationParams, p);
    model.setParameters(p);
}

void MultiConditionDataProviderHDF5::copyInputData(H5::H5File const& target)
{

    H5Ocopy(file.getId(), "/", target.getId(), "/inputData", H5P_DEFAULT, H5P_DEFAULT);
    H5Fflush(target.getId(), H5F_SCOPE_LOCAL);
}

hid_t MultiConditionDataProviderHDF5::getHdf5FileId() const { return file.getId(); }


//void MultiConditionDataProvider::printInfo() const {
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



void MultiConditionDataProviderHDF5::checkDataIntegrity() const {
    int numConditions = getNumberOfSimulationConditions();

    auto model = getModel();

    int d1, d2, d3;

    auto lock = hdf5MutexGetLock();

    assert(H5Lexists(file.getId(), hdf5MeasurementPath.c_str(), H5P_DEFAULT));
    assert(H5Lexists(file.getId(), hdf5MeasurementSigmaPath.c_str(), H5P_DEFAULT));

//    parpe::hdf5GetDatasetDimensions(file.getId(), hdf5MeasurementPath.c_str(),
//                                    3, &d1, &d2, &d3);
//    RELEASE_ASSERT(d1 >= numConditions, "");
//    RELEASE_ASSERT(d2 == model->nytrue, "");
//    RELEASE_ASSERT(d3 >= model->nt(), "");

//    parpe::hdf5GetDatasetDimensions(file.getId(), hdf5MeasurementSigmaPath.c_str(),
//                                    3, &d1, &d2, &d3);
//    RELEASE_ASSERT(d1 >= numConditions, "");
//    RELEASE_ASSERT(d2 == model->nytrue, "");
//    RELEASE_ASSERT(d3 >= model->nt(), "");

    parpe::hdf5GetDatasetDimensions(file.getId(), hdf5ConditionPath.c_str(),
                                    2, &d1, &d2);
    RELEASE_ASSERT(d1 == model->nk(), "");
}


MultiConditionDataProviderDefault::MultiConditionDataProviderDefault(std::unique_ptr<amici::Model> model, std::unique_ptr<amici::Solver> solver)
    :model(std::move(model)), solver(std::move(solver))
{

}

int MultiConditionDataProviderDefault::getNumberOfSimulationConditions() const
{
    return edata.size();
}

std::vector<int> MultiConditionDataProviderDefault::getSimulationToOptimizationParameterMapping(int  /*conditionIdx*/) const
{
    std::vector<int> mapping(model->np());
    std::iota(mapping.begin(), mapping.end(), 0);
    return mapping;
}

void MultiConditionDataProviderDefault::mapSimulationToOptimizationVariablesAddMultiply(int conditionIdx, gsl::span<double const> simulation, gsl::span<double> optimization, double coefficient) const
{
    // TODO redundant
    auto mapping = getSimulationToOptimizationParameterMapping(conditionIdx);

    for(int i = 0; i < model->np(); ++i) {
        optimization[mapping[i]] = coefficient * simulation[i];
    }
}

void MultiConditionDataProviderDefault::mapAndSetOptimizationToSimulationVariables(int conditionIdx, gsl::span<const double> optimization, gsl::span<double> simulation) const
{
    // TODO redundant
    auto mapping = getSimulationToOptimizationParameterMapping(conditionIdx);

    for(int i = 0; i < model->np(); ++i) {
        simulation[i] = optimization[mapping[i]];
    }

}

amici::ParameterScaling MultiConditionDataProviderDefault::getParameterScale(int optimizationParameterIndex) const
{
    // TODO assumes no extra optimization parameters
    return model->getParameterScale()[optimizationParameterIndex];
}


void MultiConditionDataProviderDefault::updateSimulationParameters(int  /*conditionIndex*/, gsl::span<const double> optimizationParams, amici::Model &model) const
{
    logmessage(LOGLVL_WARNING, "MultiConditionDataProviderDefault::updateSimulationParameters: No proper mapping implemented. Ensure this is correct.");
    model.setParameters(std::vector<double>(optimizationParams.begin(), optimizationParams.end()));
}

std::unique_ptr<amici::ExpData> MultiConditionDataProviderDefault::getExperimentalDataForCondition(int conditionIdx) const
{
    return std::make_unique<amici::ExpData>(edata[conditionIdx]);
}

std::vector<std::vector<double> > MultiConditionDataProviderDefault::getAllMeasurements() const
{
    std::vector<std::vector<double> > measurements;
    for(const auto& e: edata) {
        measurements.push_back(e.getObservedData());
    }
    return measurements;
}

std::vector<std::vector<double> > MultiConditionDataProviderDefault::getAllSigmas() const
{
    std::vector<std::vector<double> > sigmas;
    for(const auto& e: edata) {
        sigmas.push_back(e.getObservedDataStdDev());
    }
    return sigmas;
}


int MultiConditionDataProviderDefault::getNumOptimizationParameters() const
{
    return model->np();
}

std::unique_ptr<amici::Model> MultiConditionDataProviderDefault::getModel() const
{
    return std::unique_ptr<amici::Model>(model->clone());
}

std::unique_ptr<amici::Solver> MultiConditionDataProviderDefault::getSolver() const
{
    return std::unique_ptr<amici::Solver>(solver->clone());
}


} // namespace parpe
