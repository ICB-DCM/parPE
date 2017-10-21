#include "MultiConditionDataProvider.h"
#include "logging.h"
#include "misc.h"
#include <amici_hdf5.h>
#include <amici_misc.h>
#include <amici_model.h>
#include <cassert>
#include <cstring>
#include <edata.h>
#include <udata.h>

namespace parpe {

/**
 * @brief
 * @param hdf5Filename Filename from where to read data
 */

MultiConditionDataProvider::MultiConditionDataProvider(amici::Model *model,
                                                       std::string hdf5Filename)
    : MultiConditionDataProvider(model, hdf5Filename, "") {}

MultiConditionDataProvider::MultiConditionDataProvider(amici::Model *model,
                                                       std::string hdf5Filename,
                                                       std::string rootPath)
    : model(model) {
    hdf5LockMutex();

    H5_SAVE_ERROR_HANDLER;
    fileId = H5Fopen(hdf5Filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (fileId < 0) {
        logmessage(LOGLVL_CRITICAL,
                   "initDataProvider failed to open HDF5 file '%s'.",
                   hdf5Filename.c_str());
        printBacktrace(20);
        H5Ewalk2(H5E_DEFAULT, H5E_WALK_DOWNWARD, hdf5ErrorStackWalker_cb, NULL);
    }
    H5_RESTORE_ERROR_HANDLER;

    hdf5UnlockMutex();

    hdf5MeasurementPath = rootPath + "/measurements/y";
    hdf5MeasurementSigmaPath = rootPath + "/measurements/ysigma";
    hdf5ConditionPath = rootPath + "/fixedParameters/k";
    hdf5AmiciOptionPath = rootPath + "/amiciOptions";
    hdf5ParameterPath = rootPath + "/parameters";
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

    hdf5LockMutex();

    int d1, d2, d3;
    hdf5GetDatasetDimensions3D(fileId, hdf5MeasurementPath.c_str(), &d1, &d2, &d3);

    hdf5UnlockMutex();

    assert(d1 >= 0);

    return d1;
}

int MultiConditionDataProvider::getNumConditionSpecificParametersPerSimulation()
    const {
    hdf5LockMutex();

    int num = 0;
    int status =
        amici::AMI_HDF5_getIntScalarAttribute(fileId, hdf5ParameterPath.c_str(),
                                       "numConditionSpecificParameters", &num);
    assert(status >= 0);
    hdf5UnlockMutex();

    return num;
}

/**
 * @brief Update the contstants in AMICI UserData object. Reads a slab for the
 * given experimental conditions
 * from fixed parameters matrix.
 * @param conditionIdx Index of the experimental condition for which the
 * parameters should be taken.
 * @param udata The object to be updated.
 * @return
 */
int MultiConditionDataProvider::updateFixedSimulationParameters(
    int conditionIdx, amici::UserData &udata) const {
    hdf5LockMutex();

    H5_SAVE_ERROR_HANDLER;

    hdf5Read2DDoubleHyperslab(fileId, hdf5ConditionPath.c_str(), model->nk, 1,
                              0, conditionIdx, udata.k);

    if (H5Eget_num(H5E_DEFAULT)) {
        logmessage(LOGLVL_CRITICAL,
                   "Problem in readFixedParameters (row %d, nk %d)\n",
                   conditionIdx, model->nk);
        printBacktrace(20);
        H5Ewalk2(H5E_DEFAULT, H5E_WALK_DOWNWARD, hdf5ErrorStackWalker_cb, NULL);
        abort();
    }

    H5_RESTORE_ERROR_HANDLER;

    hdf5UnlockMutex();

    return H5Eget_num(H5E_DEFAULT);
}

std::unique_ptr<amici::ExpData> MultiConditionDataProvider::getExperimentalDataForCondition(
    int conditionIdx, const amici::UserData *udata) const {
    hdf5LockMutex();

    auto edata = std::make_unique<amici::ExpData>(udata, model);
    assert(edata && "Failed getting experimental data. Check data file.");

    hdf5Read3DDoubleHyperslab(fileId, hdf5MeasurementPath.c_str(),
                              1, edata->nytrue, edata->nt,
                              conditionIdx, 0, 0, edata->my);
    hdf5Read3DDoubleHyperslab(fileId, hdf5MeasurementSigmaPath.c_str(),
                              1, edata->nytrue, edata->nt,
                              conditionIdx, 0, 0, edata->sigmay);
    hdf5UnlockMutex();

    return edata;
}

void MultiConditionDataProvider::getOptimizationParametersLowerBounds(
    double *buffer) const {
    // TODO to HDF5
    amici::fillArray(buffer, getNumOptimizationParameters(), -2);
}

void MultiConditionDataProvider::getOptimizationParametersUpperBounds(
    double *buffer) const {
    // TODO to HDF5
    amici::fillArray(buffer, getNumOptimizationParameters(), 2);
}

int MultiConditionDataProvider::getNumOptimizationParameters() const {
    return getNumCommonParameters() +
           getNumberOfConditions() *
               getNumConditionSpecificParametersPerSimulation();
}

int MultiConditionDataProvider::getNumCommonParameters() const {
    return model->np - getNumConditionSpecificParametersPerSimulation();
}

amici::Model *MultiConditionDataProvider::getModel() const { return model; }

std::unique_ptr<amici::UserData> MultiConditionDataProvider::getUserData() const {
    // TODO: separate class for udata?
    hdf5LockMutex();

    const char *optionsObject = hdf5AmiciOptionPath.c_str();

    H5_SAVE_ERROR_HANDLER;
    auto udata = std::unique_ptr<amici::UserData>(
                amici::AMI_HDF5_readSimulationUserDataFromFileObject(
                    fileId, optionsObject, model));

    if (H5Eget_num(H5E_DEFAULT)) {
        error("Problem with reading udata\n");
        printBacktrace(20);
        H5Ewalk2(H5E_DEFAULT, H5E_WALK_DOWNWARD, hdf5ErrorStackWalker_cb, NULL);
    }
    H5_RESTORE_ERROR_HANDLER;

    hdf5UnlockMutex();

    if (udata == NULL)
        return NULL;

    // ensure parameter indices are within range //TODO: to ami_hdf
    for (int i = 0; i < udata->nplist; ++i) {
        assert(udata->plist[i] < model->np);
        assert(udata->plist[i] >= 0);
    }

    if (udata->p == NULL) {
        udata->p = new double[model->np];
    }

    if (udata->k == NULL) {
        udata->k = new double[model->nk];
    }

    // x0data is not written to HDF5 with proper dimensions
    if (udata->x0data)
        delete[] udata->x0data;
    udata->x0data = new double[model->nx]();

    return (udata);
}

std::unique_ptr<amici::UserData> MultiConditionDataProvider::getUserDataForCondition(int conditionIdx) const {
    auto udata = getUserData();

    updateFixedSimulationParameters(conditionIdx, *udata);

    return udata;
}

int MultiConditionDataProvider::
    getIndexOfFirstConditionSpecificOptimizationParameter(
        int conditionIdx) const {
    return getNumCommonParameters() +
            conditionIdx * getNumConditionSpecificParametersPerSimulation();
}

void MultiConditionDataProvider::updateSimulationParameters(int conditionIndex, const double *optimizationParams, amici::UserData *udata) const
{
    // copy all common parameters + conditionspecific parameters for first conditions to UserDaata
    memcpy(udata->p, optimizationParams, sizeof(double) * getNumCommonParameters());

    updateConditionSpecificSimulationParameters(conditionIndex, optimizationParams, udata);
}

void MultiConditionDataProvider::updateConditionSpecificSimulationParameters(
    int conditionIndex, const double *optimizationParams,
    amici::UserData *udata) const {
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

    // beginning of condition specific simulation parameters within simulation
    // parameters
    double *pConditionSpecificSimulation = &(udata->p[numCommonParams]);

    memcpy(pConditionSpecificSimulation, pConditionSpecificOptimization,
           numSpecificParams * sizeof(double));
}

hid_t MultiConditionDataProvider::getHdf5FileId() const { return fileId; }

MultiConditionDataProvider::~MultiConditionDataProvider() {
    if (fileId > 0) {
        H5_SAVE_ERROR_HANDLER;
        herr_t status = H5Fclose(fileId);

        if (status < 0) {
            error("closeDataProvider failed to close HDF5 file.");
            printBacktrace(20);
            H5Ewalk2(H5E_DEFAULT, H5E_WALK_DOWNWARD, hdf5ErrorStackWalker_cb,
                     NULL);
        }
        H5_RESTORE_ERROR_HANDLER;
    }
}

void JobIdentifier::print() const {
    printf("%d.%d.%d.%d", idxMultiStart, idxLocalOptimization,
           idxLocalOptimizationIteration, idxConditions);
}

void JobIdentifier::sprint(char *buffer) const {
    sprintf(buffer, "%d.%d.%d.%d", idxMultiStart, idxLocalOptimization,
            idxLocalOptimizationIteration, idxConditions);
}

} // namespace parpe
