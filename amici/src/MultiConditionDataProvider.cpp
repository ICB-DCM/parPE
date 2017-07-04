#include "MultiConditionDataProvider.h"
#include "logging.h"
#include "misc.h"
#include <assert.h>
#include "ami_hdf5.h"

UserData getUserData();
// alias because getUserData is shadowed in MultiConditionDataProvider
inline UserData getModelUserData() { return getUserData(); }

void printJobIdentifierifier(JobIdentifier id)
{
    printf("%d.%d.%d.%d", id.idxMultiStart, id.idxLocalOptimization, id.idxLocalOptimizationIteration, id.idxConditions);
}

void sprintJobIdentifier(char *buffer, JobIdentifier id)
{
    sprintf(buffer, "%d.%d.%d.%d", id.idxMultiStart, id.idxLocalOptimization, id.idxLocalOptimizationIteration, id.idxConditions);
}


/**
 * @brief
 * @param hdf5Filename Filename from where to read data
 */

MultiConditionDataProvider::MultiConditionDataProvider(const char *hdf5Filename) : modelDims(getModelUserData())
{
    H5_SAVE_ERROR_HANDLER;
    fileId = H5Fopen(hdf5Filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    if(fileId < 0) {
        logmessage(LOGLVL_CRITICAL, "initDataProvider failed to open HDF5 file '%s'.", hdf5Filename);
        printBacktrace(20);
        H5Ewalk2(H5E_DEFAULT, H5E_WALK_DOWNWARD, hdf5ErrorStackWalker_cb, NULL);
    }
    H5_RESTORE_ERROR_HANDLER;

    hdf5MeasurementPath = "/measurements/y";
    hdf5MeasurementSigmaPath = "/measurements/ysigma";
    hdf5ConditionPath = "/fixedParameters/k";
}


/**
 * @brief Get the number of simulations required for objective function evaluation. Currently, this amounts to the number
 * of conditions present in the data.
 * @return
 */
int MultiConditionDataProvider::getNumberOfConditions()
{
    // TODO: add additional layer for selecten of condition indices (for testing and later for minibatch)
    // -> won't need different file for testing/validation splits
    // TODO: cache

    hid_t dset = H5Dopen(fileId, hdf5MeasurementPath.c_str(), H5P_DEFAULT);
    hid_t dspace = H5Dget_space(dset);
    const int ndims = H5Sget_simple_extent_ndims(dspace);
    assert(ndims == 2); // nObservables * nConditions (TODO *nt)
    hsize_t dims[ndims];
    H5Sget_simple_extent_dims(dspace, dims, NULL);

    return dims[1];
}

int MultiConditionDataProvider::getNumConditionSpecificParametersPerSimulation()
{
    return AMI_HDF5_getIntScalarAttribute(fileId, "/parameters", "numConditionSpecificParameters");
}


/**
 * @brief Update the contstants in AMICI UserData object. Reads a slab for the given experimental conditions
 * from fixed parameters matrix.
 * @param conditionIdx Index of the experimental condition for which the parameters should be taken.
 * @param udata The object to be updated.
 * @return
 */
int MultiConditionDataProvider::updateFixedSimulationParameters(int conditionIdx, UserData *udata)
{
    hdf5LockMutex();

    H5_SAVE_ERROR_HANDLER;

    hdf5Read2DDoubleHyperslab(fileId, hdf5ConditionPath.c_str(), udata->nk, 1, 0, conditionIdx, udata->k);

    if(H5Eget_num(H5E_DEFAULT)) {
        logmessage(LOGLVL_CRITICAL, "Problem in readFixedParameters (row %d, nk %d)\n", conditionIdx, udata->nk);
        printBacktrace(20);
        H5Ewalk2(H5E_DEFAULT, H5E_WALK_DOWNWARD, hdf5ErrorStackWalker_cb, NULL);
        abort();
    }

    H5_RESTORE_ERROR_HANDLER;

    hdf5UnlockMutex();

    return H5Eget_num(H5E_DEFAULT);

}


/**
 * @brief Reads data required for simulation of a specific experimental condition. Creates ExpData and updates UserData object.
 * @param dpath
 * @param udata
 * @return
 */

ExpData *MultiConditionDataProvider::getExperimentalDataForExperimentAndUpdateUserData(int conditionIdx, UserData *udata)
{

    updateFixedSimulationParameters(conditionIdx, udata);


    return getExperimentalDataForCondition(conditionIdx);

}

ExpData *MultiConditionDataProvider::getExperimentalDataForCondition(int conditionIdx)
{
    hdf5LockMutex();

    ExpData *edata = new ExpData(&modelDims);

    hdf5Read2DDoubleHyperslab(fileId, hdf5MeasurementPath.c_str(),      modelDims.ny, 1, 0, conditionIdx, edata->my);
    hdf5Read2DDoubleHyperslab(fileId, hdf5MeasurementSigmaPath.c_str(), modelDims.ny, 1, 0, conditionIdx, edata->sigmay);

    hdf5UnlockMutex();

    return edata;

}


void MultiConditionDataProvider::getOptimizationParametersLowerBounds(double *buffer)
{
    // TODO to HDF5
    fillArray(buffer, modelDims.np, -2);
}


void MultiConditionDataProvider::getOptimizationParametersUpperBounds(double *buffer)
{
    // TODO to HDF5
    fillArray(buffer, modelDims.np, 2);

}


int MultiConditionDataProvider::getNumOptimizationParameters()
{
    return getNumCommonParameters() + getNumberOfConditions() * getNumConditionSpecificParametersPerSimulation();

}


int MultiConditionDataProvider::getNumCommonParameters()
{
    return modelDims.np - getNumConditionSpecificParametersPerSimulation();

}

UserData MultiConditionDataProvider::getModelDims()
{
    return modelDims;
}

UserData *MultiConditionDataProvider::getUserData()
{
    // TODO: separate class for udata?
    hdf5LockMutex();

    const char* optionsObject = "/amiciOptions";

    H5_SAVE_ERROR_HANDLER;
    UserData *udata = AMI_HDF5_readSimulationUserDataFromFileObject(fileId, optionsObject);

    if(H5Eget_num(H5E_DEFAULT)) {
        error("Problem with reading udata\n");
        printBacktrace(20);
        H5Ewalk2(H5E_DEFAULT, H5E_WALK_DOWNWARD, hdf5ErrorStackWalker_cb, NULL);
    }
    H5_RESTORE_ERROR_HANDLER;

    hdf5UnlockMutex();

    if(udata == NULL)
        return NULL;

    // ensure parameter indices are within range //TODO: to ami_hdf
    for(int i = 0; i < udata->nplist; ++i) {
        assert(udata->plist[i] < udata->np);
        assert(udata->plist[i] >= 0);
    }

    if(udata->p == NULL) {
        udata->p = new double[udata->np];
    }

    if(udata->k == NULL) {
        udata->k = new double[udata->nk];
    }

    // x0data is not written to HDF5 with proper dimensions
    if(udata->x0data)
        delete[] udata->x0data;
    udata->x0data = new double[udata->nx]();

    udata->pscale = AMICI_SCALING_LOG10;

    return(udata);
}

UserData *MultiConditionDataProvider::getUserDataForCondition(int conditionIdx)
{
    UserData *udata = getUserData();
    updateFixedSimulationParameters(conditionIdx, udata);
    return udata;
}

int MultiConditionDataProvider::getIndexOfFirstConditionSpecificOptimizationParameter(int conditionIdx)
{
    return getNumCommonParameters() + conditionIdx * getNumConditionSpecificParametersPerSimulation();
}


MultiConditionDataProvider::~MultiConditionDataProvider()
{
    H5_SAVE_ERROR_HANDLER;
    herr_t status = H5Fclose(fileId);

    if(status< 0) {
        error("closeDataProvider failed to close HDF5 file.");
        printBacktrace(20);
        H5Ewalk2(H5E_DEFAULT, H5E_WALK_DOWNWARD, hdf5ErrorStackWalker_cb, NULL);
    }
    H5_RESTORE_ERROR_HANDLER;
}
