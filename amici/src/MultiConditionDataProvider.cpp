#include "MultiConditionDataProvider.h"
#include "logging.h"
#include "misc.h"
#include <amici_hdf5.h>
#include <amici_misc.h>
#include <amici_model.h>
#include <cassert>
#include <cstring>
#include <edata.h>

namespace parpe {

/**
 * @brief
 * @param hdf5Filename Filename from where to read data
 */

MultiConditionDataProvider::MultiConditionDataProvider(std::unique_ptr<amici::Model> model,
                                                       std::string hdf5Filename)
    : MultiConditionDataProvider(std::move(model), hdf5Filename, "") {}

MultiConditionDataProvider::MultiConditionDataProvider(std::unique_ptr<amici::Model> model,
                                                       std::string hdf5Filename,
                                                       std::string rootPath)
    : model(std::move(model)), rootPath(rootPath) {

    auto lock = hdf5MutexGetLock();

    H5_SAVE_ERROR_HANDLER;
    try {
        file = H5::H5File(hdf5Filename.c_str(), H5F_ACC_RDONLY);
        fileId = file.getId();
    } catch (...) {
        logmessage(LOGLVL_CRITICAL,
                   "initDataProvider failed to open HDF5 file '%s'.",
                   hdf5Filename.c_str());
        printBacktrace(20);
        H5Ewalk2(H5E_DEFAULT, H5E_WALK_DOWNWARD, hdf5ErrorStackWalker_cb, NULL);
        throw(HDF5Exception());
    }
    H5_RESTORE_ERROR_HANDLER;

    hdf5MeasurementPath = rootPath + "/measurements/y";
    hdf5MeasurementSigmaPath = rootPath + "/measurements/ysigma";
    hdf5ConditionPath = rootPath + "/fixedParameters/k";
    hdf5AmiciOptionPath = rootPath + "/amiciOptions";
    hdf5ParameterPath = rootPath + "/parameters";
    hdf5ParameterMinPath = hdf5ParameterPath + "/lowerBound";
    hdf5ParameterMaxPath = hdf5ParameterPath + "/upperBound";

    amici::readModelDataFromHDF5(fileId, *this->model, hdf5AmiciOptionPath.c_str());
}

/**
 * @brief Get the number of simulations required for objective function
 * evaluation. Currently, this amounts to the number
 * of conditions present in the data.
 * @return
 */
int MultiConditionDataProvider::getNumberOfConditions() const {
    // TODO: add additional layer for selecten of condition indices (for testing
    // and later for minibatch)
    // -> won't need different file for testing/validation splits
    // TODO: cache

    auto lock = hdf5MutexGetLock();

    int d1, d2, d3;
    hdf5GetDatasetDimensions(fileId, hdf5MeasurementPath.c_str(), 3, &d1, &d2, &d3);

    assert(d1 >= 0);

    return d1;
}

int MultiConditionDataProvider::getNumConditionSpecificParametersPerSimulation()
    const {
    auto lock = hdf5MutexGetLock();

    int num = 0;
    int status =
        amici::AMI_HDF5_getIntScalarAttribute(fileId, hdf5ParameterPath.c_str(),
                                       "numConditionSpecificParameters", &num);
    assert(status >= 0);

    return num;
}

/**
 * @brief Update the contstants in AMICI UserData object. Reads a slab for the
 * given experimental conditions
 * from fixed parameters matrix.
 * @param conditionIdx Index of the experimental condition for which the
 * parameters should be taken.
 * @param udata The object to be updated.
 * @return On success zero, non-zero on failure
 */
int MultiConditionDataProvider::updateFixedSimulationParameters(int conditionIdx, amici::Model &model) const {
    auto lock = hdf5MutexGetLock();

    H5_SAVE_ERROR_HANDLER;

    std::vector<double> buf(model.nk());
    hdf5Read2DDoubleHyperslab(fileId, hdf5ConditionPath.c_str(), model.nk(), 1,
                              0, conditionIdx, buf.data());
    model.setFixedParameters(buf);

    if (H5Eget_num(H5E_DEFAULT)) {
        logmessage(LOGLVL_CRITICAL,
                   "Problem in readFixedParameters (row %d, nk %d)\n",
                   conditionIdx, model.nk());
        printBacktrace(20);
        H5Ewalk2(H5E_DEFAULT, H5E_WALK_DOWNWARD, hdf5ErrorStackWalker_cb, NULL);
        abort();
    }

    H5_RESTORE_ERROR_HANDLER;

    return H5Eget_num(H5E_DEFAULT);
}

std::unique_ptr<amici::ExpData> MultiConditionDataProvider::getExperimentalDataForCondition(
    int conditionIdx) const {
    auto lock = hdf5MutexGetLock();

    auto edata = std::make_unique<amici::ExpData>(model.get());
    assert(edata && "Failed getting experimental data. Check data file.");

    hdf5Read3DDoubleHyperslab(fileId, hdf5MeasurementPath.c_str(),
                              1, edata->nytrue, edata->nt,
                              conditionIdx, 0, 0, edata->my.data());
    hdf5Read3DDoubleHyperslab(fileId, hdf5MeasurementSigmaPath.c_str(),
                              1, edata->nytrue, edata->nt,
                              conditionIdx, 0, 0, edata->sigmay.data());

    return edata;
}

void MultiConditionDataProvider::getOptimizationParametersLowerBounds(
    double *buffer) const {

    auto dataset = file.openDataSet(hdf5ParameterMinPath);

    auto dataspace = dataset.getSpace();
    RELEASE_ASSERT(dataspace.getSimpleExtentNdims() == 1, "hdf5ParameterMinPath dimensions dont match");
    hsize_t dim = 0;
    dataspace.getSimpleExtentDims(&dim);
    RELEASE_ASSERT(dim == (unsigned) getNumOptimizationParameters(), "hdf5ParameterMinPath dimensions dont match");

    dataset.read(buffer, H5::PredType::NATIVE_DOUBLE);
}

void MultiConditionDataProvider::getOptimizationParametersUpperBounds(
    double *buffer) const {
    auto dataset = file.openDataSet(hdf5ParameterMaxPath);

    auto dataspace = dataset.getSpace();
    RELEASE_ASSERT(dataspace.getSimpleExtentNdims() == 1, "hdf5ParameterMaxPath dimensions dont match");
    hsize_t dim = 0;
    dataspace.getSimpleExtentDims(&dim);
    RELEASE_ASSERT(dim == (unsigned) getNumOptimizationParameters(), "hdf5ParameterMaxPath dimensions dont match");

    dataset.read(buffer, H5::PredType::NATIVE_DOUBLE);
}

int MultiConditionDataProvider::getNumOptimizationParameters() const {
    return getNumCommonParameters() +
           getNumberOfConditions() *
               getNumConditionSpecificParametersPerSimulation();
}

int MultiConditionDataProvider::getNumCommonParameters() const {
    return model->np() - getNumConditionSpecificParametersPerSimulation();
}

std::unique_ptr<amici::Model> MultiConditionDataProvider::getModel() const { return std::unique_ptr<amici::Model>(model->clone()); }

std::unique_ptr<amici::Solver> MultiConditionDataProvider::getSolver() const
{
    auto solver = model->getSolver();
    amici::readSolverSettingsFromHDF5(fileId, *solver, hdf5AmiciOptionPath.c_str());
    return solver;
}


std::unique_ptr<amici::Model> MultiConditionDataProvider::getModelForCondition(int conditionIdx) const {
    // TODO: RENAME. obsolete?
    auto newModel = getModel();

    updateFixedSimulationParameters(conditionIdx, *newModel);

    return newModel;
}

int MultiConditionDataProvider::
    getIndexOfFirstConditionSpecificOptimizationParameter(
        int conditionIdx) const {
    return getNumCommonParameters() +
            conditionIdx * getNumConditionSpecificParametersPerSimulation();
}

void MultiConditionDataProvider::updateSimulationParameters(int conditionIndex, const double *optimizationParams, amici::Model &model) const
{
    // copy all common parameters + conditionspecific parameters for first conditions to UserDaata
    auto p = model.getParameters();
    std::copy(optimizationParams, optimizationParams + getNumCommonParameters(), p.data());
    model.setParameters(p);
    updateConditionSpecificSimulationParameters(conditionIndex, optimizationParams, model);
}

void MultiConditionDataProvider::updateConditionSpecificSimulationParameters(int conditionIndex, const double *optimizationParams,
    amici::Model& model) const {
    /* Optimization parameters are [commonParameters,
     * condition1SpecificParameters, condition2SpecificParameters, ...]
     * number of condition specific parameters is the same for all cell lines.
     * Simulation parameters are [commonParameters,
     * currentCelllineSpecificParameter]
     */


    const int numCommonParams = getNumCommonParameters();
    const int numSpecificParams =
        getNumConditionSpecificParametersPerSimulation();

    // beginning of condition specific simulation parameters within optimization
    // parameters
    const double *pConditionSpecificOptimization =
        &optimizationParams
            [getIndexOfFirstConditionSpecificOptimizationParameter(
                conditionIndex)];

    auto p = model.getParameters();
    // beginning of condition specific simulation parameters within simulation
    // parameters
    double *pConditionSpecificSimulation = &(p.data()[numCommonParams]);

    std::copy(pConditionSpecificOptimization, pConditionSpecificOptimization + numSpecificParams, pConditionSpecificSimulation);
    model.setParameters(p);

}

void MultiConditionDataProvider::copyInputData(H5::H5File target)
{
    H5Ocopy(fileId, "/", target.getId(), "/inputData", H5P_DEFAULT, H5P_DEFAULT);
    H5Fflush(target.getId(), H5F_SCOPE_LOCAL);
}

hid_t MultiConditionDataProvider::getHdf5FileId() const { return fileId; }


void JobIdentifier::print() const {
    printf("%d.%d.%d.%d", idxMultiStart, idxLocalOptimization,
           idxLocalOptimizationIteration, idxConditions);
}

void JobIdentifier::sprint(char *buffer) const {
    sprintf(buffer, "%d.%d.%d.%d", idxMultiStart, idxLocalOptimization,
            idxLocalOptimizationIteration, idxConditions);
}

} // namespace parpe
