#include "resultwriter.h"
#include "misc.h"
#include <string.h>
#include <assert.h>
#include <pthread.h>

#define H5_SAVE_ERROR_HANDLER   herr_t (*old_func)(void*); \
                                void *old_client_data; \
                                H5Eget_auto1(&old_func, &old_client_data); \
                                H5Eset_auto1(NULL, NULL);
#define H5_RESTORE_ERROR_HANDLER H5Eset_auto1(old_func, old_client_data);

static hid_t file_id = 0;
extern pthread_mutex_t mutexHDF;


// TODO
static bool hdf5DatasetExists(hid_t file_id, const char *datasetName);

static bool hdf5GroupExists(hid_t file_id, const char *groupName);

static void hdf5CreateGroup(hid_t file_id, const char *groupPath, bool recursively);

static void hdf5CreateExtendableDouble2DArray(hid_t file_id, const char* datasetPath, int stride);

static void hdf5CreateExtendableInt2DArray(hid_t file_id, const char* datasetPath, int stride);

static void hdf5CreateExtendableDouble3DArray(hid_t file_id, const char *datasetPath, int stride1, int stride2);

static void hdf5Extend2ndDimensionAndWriteToDouble2DArray(hid_t file_id, const char *datasetPath, const double *buffer);

static void hdf5Extend2ndDimensionAndWriteToInt2DArray(hid_t file_id, const char *datasetPath, const int *buffer);

static void hdf5CreateOrExtendAndWriteToDouble2DArray(hid_t file_id, const char *parentPath, const char *datasetName, const double *buffer, int stride);

static void hdf5CreateOrExtendAndWriteToInt2DArray(hid_t file_id, const char *parentPath, const char *datasetName, const int *buffer, int stride);

static void hdf5CreateOrExtendAndWriteToDouble3DArray(hid_t file_id, const char *parentPath, const char *datasetName, const double *buffer, int stride1, int stride2);

static char *myStringCat(const char *first, const char *second);


int initResultHDFFile(const char *filename)
{
    H5_SAVE_ERROR_HANDLER;

    file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT); // overwrites TODO: backup
    if(file_id < 0) {
        logmessage(LOGLVL_CRITICAL, "Cannot open result file '%s'", filename);
        H5Eprint(H5E_DEFAULT, stderr);
    } else {
        logmessage(LOGLVL_INFO, "Writing results to '%s'", filename);
    }

    H5_RESTORE_ERROR_HANDLER;

    return file_id < 0;
}


void closeResultHDFFile()
{
    H5_SAVE_ERROR_HANDLER
    herr_t status = H5Fclose(file_id);

    if(status< 0) {
        error("closeResultHDFFile failed to close HDF5 file.");
        H5Eprint(H5E_DEFAULT, stderr);
    }
    H5_RESTORE_ERROR_HANDLER
}

void flushResultWriter()
{
    pthread_mutex_lock(&mutexHDF);
    H5Fflush(file_id, H5F_SCOPE_LOCAL);
    pthread_mutex_unlock(&mutexHDF);
}

bool hdf5DatasetExists(hid_t file_id, const char *datasetName)
{
    pthread_mutex_lock(&mutexHDF);
    bool exists = H5Lexists(file_id, datasetName, H5P_DEFAULT) > 0;
    pthread_mutex_unlock(&mutexHDF);

    return exists;
}

bool hdf5GroupExists(hid_t file_id, const char *groupName)
{
    pthread_mutex_lock(&mutexHDF);

    // switch off error handler, check existance and reenable
    herr_t (*old_func)(void*);
    void *old_client_data;
    H5Eget_auto1(&old_func, &old_client_data);
    H5Eset_auto1(NULL, NULL);

    herr_t status = H5Gget_objinfo(file_id, groupName, 0, NULL);

    H5Eset_auto1(old_func, old_client_data);

    pthread_mutex_unlock(&mutexHDF);

    return status >= 0;
}

void hdf5CreateGroup(hid_t file_id, const char *groupPath, bool recursively)
{
    hid_t groupCreationPropertyList = H5P_DEFAULT;

    pthread_mutex_lock(&mutexHDF);

    if(recursively) {
        groupCreationPropertyList = H5Pcreate(H5P_LINK_CREATE);
        H5Pset_create_intermediate_group (groupCreationPropertyList, 1);
    }

    hid_t group = H5Gcreate(file_id, groupPath, groupCreationPropertyList, H5P_DEFAULT, H5P_DEFAULT);
    if(group < 0)
        error("Failed to create group in hdf5CreateGroup");
    H5Gclose(group);

    pthread_mutex_unlock(&mutexHDF);
}

void hdf5CreateExtendableDouble2DArray(hid_t file_id, const char *datasetPath, int stride)
{
    int rank = 2;
    hsize_t initialDimensions[2]  = {stride, 0};
    hsize_t maximumDimensions[2]  = {stride, H5S_UNLIMITED};

    pthread_mutex_lock(&mutexHDF);

    hid_t dataspace = H5Screate_simple(rank, initialDimensions, maximumDimensions);

    // need chunking for extendable dataset
    hsize_t chunkDimensions[2] = {stride, 1};
    hid_t datasetCreationProperty = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(datasetCreationProperty, rank, chunkDimensions);

    hid_t dataset = H5Dcreate2(file_id, datasetPath, H5T_NATIVE_DOUBLE, dataspace,
                               H5P_DEFAULT, datasetCreationProperty, H5P_DEFAULT);

    H5Dclose(dataset);
    H5Sclose(dataspace);

    pthread_mutex_unlock(&mutexHDF);
}

void hdf5Extend2ndDimensionAndWriteToDouble2DArray(hid_t file_id, const char *datasetPath, const double *buffer)
{
    pthread_mutex_lock(&mutexHDF);

    hid_t dataset = H5Dopen2(file_id, datasetPath, H5P_DEFAULT);

    // extend
    hid_t filespace = H5Dget_space(dataset);
    int rank = H5Sget_simple_extent_ndims(filespace);
    assert(rank == 2);

    hsize_t currentDimensions[2];
    H5Sget_simple_extent_dims(filespace, currentDimensions, NULL);

    hsize_t newDimensions[2]  = {currentDimensions[0], currentDimensions[1] + 1};
    herr_t status = H5Dset_extent(dataset, newDimensions);

    filespace = H5Dget_space(dataset);
    hsize_t offset[2] = {0, currentDimensions[1]};
    hsize_t slabsize[2] = {currentDimensions[0], 1};

    status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, slabsize, NULL);

    hid_t memspace = H5Screate_simple(rank, slabsize, NULL);

    status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, memspace, filespace, H5P_DEFAULT, buffer);

    if(status < 0)
        error("Failed to write data in hdf5Extend2ndDimensionAndWriteToDouble2DArray");

    H5Dclose(dataset);
    H5Sclose(filespace);
    H5Sclose(memspace);

    pthread_mutex_unlock(&mutexHDF);
}

void hdf5Extend3rdDimensionAndWriteToDouble3DArray(hid_t file_id, const char *datasetPath, const double *buffer)
{
    pthread_mutex_lock(&mutexHDF);

    hid_t dataset = H5Dopen2(file_id, datasetPath, H5P_DEFAULT);

    // extend
    hid_t filespace = H5Dget_space(dataset);
    int rank = H5Sget_simple_extent_ndims(filespace);
    assert(rank == 3);

    hsize_t currentDimensions[3];
    H5Sget_simple_extent_dims(filespace, currentDimensions, NULL);

    hsize_t newDimensions[3]  = {currentDimensions[0], currentDimensions[1], currentDimensions[2] + 1};
    herr_t status = H5Dset_extent(dataset, newDimensions);

    filespace = H5Dget_space(dataset);
    hsize_t offset[3] = {0, 0, currentDimensions[2]};
    hsize_t slabsize[3] = {currentDimensions[0], currentDimensions[1], 1};

    status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, slabsize, NULL);

    hid_t memspace = H5Screate_simple(rank, slabsize, NULL);

    status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, memspace, filespace, H5P_DEFAULT, buffer);

    if(status < 0)
        error("Failed to write data in hdf5Extend3rdDimensionAndWriteToDouble3DArray");

    H5Dclose(dataset);
    H5Sclose(filespace);
    H5Sclose(memspace);

    pthread_mutex_unlock(&mutexHDF);
}


char *myStringCat(const char *first, const char *second)
{
    char *concatenation = malloc(sizeof(char) * (strlen(first) + strlen(second) + 1));
    strcpy(concatenation, first);
    strcat(concatenation, second);
    return concatenation;
}

void logLocalOptimizerObjectiveFunctionEvaluation(datapath path, int numFunctionCalls, double *theta, double objectiveFunctionValue, double timeElapsedInSeconds, int nTheta)
{
    char fullGroupPath[100];
    sprintf(fullGroupPath, "/crossvalidations/%d/multistarts/%d/objFunEval/", path.idxMultiStart, path.idxLocalOptimization);

    hdf5CreateOrExtendAndWriteToDouble2DArray(file_id, fullGroupPath, "negLogLikelihood", &objectiveFunctionValue, 1);
    hdf5CreateOrExtendAndWriteToDouble2DArray(file_id, fullGroupPath, "p", theta, nTheta);
    hdf5CreateOrExtendAndWriteToDouble2DArray(file_id, fullGroupPath, "evalFTime", &timeElapsedInSeconds, 1);

    flushResultWriter();
}


void logLocalOptimizerObjectiveFunctionGradientEvaluation(datapath path, int numFunctionCalls, double *theta, double objectiveFunctionValue, const double *gradient, double timeElapsedInSeconds, int nTheta)
{
    char fullGroupPath[100];
    sprintf(fullGroupPath, "/crossvalidations/%d/multistarts/%d/objFunGradEval/", path.idxMultiStart, path.idxLocalOptimization);

    hdf5CreateOrExtendAndWriteToDouble2DArray(file_id, fullGroupPath, "negLogLikelihood", &objectiveFunctionValue, 1);
    hdf5CreateOrExtendAndWriteToDouble2DArray(file_id, fullGroupPath, "negLogLikelihoodGradient", gradient, nTheta);
    hdf5CreateOrExtendAndWriteToDouble2DArray(file_id, fullGroupPath, "p", theta, nTheta);
    hdf5CreateOrExtendAndWriteToDouble2DArray(file_id, fullGroupPath, "evalFTime", &timeElapsedInSeconds, 1);

    flushResultWriter();
}

void logLocalOptimizerIteration(datapath path, int numIterations, double *theta, double objectiveFunctionValue, const double *gradient, double timeElapsedInSeconds, int nTheta, int alg_mod, double inf_pr, double inf_du, double mu, double d_norm, double regularization_size, double alpha_du, double alpha_pr, int ls_trials)
{
    char fullGroupPath[100];
    sprintf(fullGroupPath, "/crossvalidations/%d/multistarts/%d/optimizerIterations/", path.idxMultiStart, path.idxLocalOptimization);

    hdf5CreateOrExtendAndWriteToDouble2DArray(file_id, fullGroupPath, "negLogLikelihood", &objectiveFunctionValue, 1);
    hdf5CreateOrExtendAndWriteToDouble2DArray(file_id, fullGroupPath, "negLogLikelihoodGradient", gradient, nTheta);
    hdf5CreateOrExtendAndWriteToDouble2DArray(file_id, fullGroupPath, "p", theta, nTheta);
    hdf5CreateOrExtendAndWriteToDouble2DArray(file_id, fullGroupPath, "evalFTime", &timeElapsedInSeconds, 1);

    hdf5CreateOrExtendAndWriteToInt2DArray(file_id, fullGroupPath, "numIterations", &numIterations, 1);
    hdf5CreateOrExtendAndWriteToInt2DArray(file_id, fullGroupPath, "alg_mod", &alg_mod, 1);
    hdf5CreateOrExtendAndWriteToDouble2DArray(file_id, fullGroupPath, "inf_pr", &inf_pr, 1);
    hdf5CreateOrExtendAndWriteToDouble2DArray(file_id, fullGroupPath, "inf_du", &inf_du, 1);
    hdf5CreateOrExtendAndWriteToDouble2DArray(file_id, fullGroupPath, "mu", &mu, 1);
    hdf5CreateOrExtendAndWriteToDouble2DArray(file_id, fullGroupPath, "d_norm", &d_norm, 1);
    hdf5CreateOrExtendAndWriteToDouble2DArray(file_id, fullGroupPath, "regularization_size", &regularization_size, 1);
    hdf5CreateOrExtendAndWriteToDouble2DArray(file_id, fullGroupPath, "alpha_du", &alpha_du, 1);
    hdf5CreateOrExtendAndWriteToDouble2DArray(file_id, fullGroupPath, "alpha_pr", &alpha_pr, 1);
    hdf5CreateOrExtendAndWriteToInt2DArray(file_id, fullGroupPath, "ls_trials", &ls_trials, 1);

    flushResultWriter();
}

void logSimulation(datapath path, double llh, const double *gradient, double timeElapsedInSeconds, int nTheta, int numStates, double *states, double *stateSensi, double *y, int jobId, int iterationsUntilSteadystate)
{
    char fullGroupPath[200];
    // TODO doesnt have to be extendable in current implementation
    sprintf(fullGroupPath, "/crossvalidations/%d/multistarts/%d/iteration/%d/genotype/%d/experiment/%d/",
            path.idxMultiStart, path.idxLocalOptimization, path.idxLocalOptimizationIteration, path.idxGenotype, path.idxExperiment);

    hdf5CreateOrExtendAndWriteToDouble2DArray(file_id, fullGroupPath, "negLogLikelihood", &llh, 1);
    hdf5CreateOrExtendAndWriteToInt2DArray(file_id, fullGroupPath, "jobId", &jobId, 1);
    hdf5CreateOrExtendAndWriteToDouble2DArray(file_id, fullGroupPath, "negLogLikelihoodGradient", gradient, nTheta);
    //hdf5CreateOrExtendAndWriteToDouble2DArray(file_id, fullGroupPath, "p", theta, nTheta);
    hdf5CreateOrExtendAndWriteToDouble2DArray(file_id, fullGroupPath, "evalFTime", &timeElapsedInSeconds, 1);
    hdf5CreateOrExtendAndWriteToDouble2DArray(file_id, fullGroupPath, "X", states, numStates);
    hdf5CreateOrExtendAndWriteToDouble2DArray(file_id, fullGroupPath, "Y", y, 1);
    hdf5CreateOrExtendAndWriteToInt2DArray(file_id, fullGroupPath, "iterationsUntilSteadystate", &iterationsUntilSteadystate, 1);

    if(stateSensi)
        hdf5CreateOrExtendAndWriteToDouble3DArray(file_id, fullGroupPath, "sX", stateSensi, numStates, nTheta);

    flushResultWriter();
}


void hdf5CreateOrExtendAndWriteToDouble2DArray(hid_t file_id, const char *parentPath, const char *datasetName, const double *buffer, int stride)
{
    pthread_mutex_lock(&mutexHDF);

    if(!hdf5GroupExists(file_id, parentPath)) {
        hdf5CreateGroup(file_id, parentPath, true);
    }

    char *fullDatasetPath = myStringCat(parentPath, datasetName);

    if(!hdf5DatasetExists(file_id, fullDatasetPath)) {
        hdf5CreateExtendableDouble2DArray(file_id, fullDatasetPath, stride);
    }

    hdf5Extend2ndDimensionAndWriteToDouble2DArray(file_id, fullDatasetPath, buffer);

    free(fullDatasetPath);

    pthread_mutex_unlock(&mutexHDF);
}

void hdf5CreateOrExtendAndWriteToDouble3DArray(hid_t file_id, const char *parentPath, const char *datasetName, const double *buffer, int stride1, int stride2)
{
    pthread_mutex_lock(&mutexHDF);

    if(!hdf5GroupExists(file_id, parentPath)) {
        hdf5CreateGroup(file_id, parentPath, true);
    }

    char *fullDatasetPath = myStringCat(parentPath, datasetName);

    if(!hdf5DatasetExists(file_id, fullDatasetPath)) {
        hdf5CreateExtendableDouble3DArray(file_id, fullDatasetPath, stride1, stride2);
    }

    hdf5Extend3rdDimensionAndWriteToDouble3DArray(file_id, fullDatasetPath, buffer);

    pthread_mutex_unlock(&mutexHDF);

    free(fullDatasetPath);
}


void hdf5CreateOrExtendAndWriteToInt2DArray(hid_t file_id, const char *parentPath, const char *datasetName, const int *buffer, int stride)
{
    pthread_mutex_lock(&mutexHDF);

    if(!hdf5GroupExists(file_id, parentPath)) {
        hdf5CreateGroup(file_id, parentPath, true);
    }

    char *fullDatasetPath = myStringCat(parentPath, datasetName);

    if(!hdf5DatasetExists(file_id, fullDatasetPath)) {
        hdf5CreateExtendableInt2DArray(file_id, fullDatasetPath, stride);
    }

    hdf5Extend2ndDimensionAndWriteToInt2DArray(file_id, fullDatasetPath, buffer);

    pthread_mutex_unlock(&mutexHDF);

    free(fullDatasetPath);
}

void hdf5Extend2ndDimensionAndWriteToInt2DArray(hid_t file_id, const char *datasetPath, const int *buffer)
{
    pthread_mutex_lock(&mutexHDF);

    hid_t dataset = H5Dopen2(file_id, datasetPath, H5P_DEFAULT);

    // extend
    hid_t filespace = H5Dget_space(dataset);
    int rank = H5Sget_simple_extent_ndims(filespace);
    assert(rank == 2);

    hsize_t currentDimensions[2];
    H5Sget_simple_extent_dims(filespace, currentDimensions, NULL);

    hsize_t newDimensions[2]  = {currentDimensions[0], currentDimensions[1] + 1};
    herr_t status = H5Dset_extent(dataset, newDimensions);

    filespace = H5Dget_space(dataset);
    hsize_t offset[2] = {0, currentDimensions[1]};
    hsize_t slabsize[2] = {currentDimensions[0], 1};

    status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, slabsize, NULL);

    hid_t memspace = H5Screate_simple(rank, slabsize, NULL);

    status = H5Dwrite(dataset, H5T_NATIVE_INT, memspace, filespace, H5P_DEFAULT, buffer);

    if(status < 0)
        error("Error writing data in hdf5Extend2ndDimensionAndWriteToInt2DArray.");

    H5Dclose(dataset);
    H5Sclose(filespace);
    H5Sclose(memspace);

    pthread_mutex_unlock(&mutexHDF);
}

void hdf5CreateExtendableInt2DArray(hid_t file_id, const char *datasetPath, int stride)
{
    pthread_mutex_lock(&mutexHDF);

    int rank = 2;
    hsize_t initialDimensions[2]  = {stride, 0};
    hsize_t maximumDimensions[2]  = {stride, H5S_UNLIMITED};

    hid_t dataspace = H5Screate_simple(rank, initialDimensions, maximumDimensions);

    // need chunking for extendable dataset
    hsize_t chunkDimensions[2] = {stride, 1};
    hid_t datasetCreationProperty = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(datasetCreationProperty, rank, chunkDimensions);

    hid_t dataset = H5Dcreate2(file_id, datasetPath, H5T_NATIVE_INT, dataspace,
                               H5P_DEFAULT, datasetCreationProperty, H5P_DEFAULT);

    H5Dclose(dataset);
    H5Sclose(dataspace);

    pthread_mutex_unlock(&mutexHDF);
}


void hdf5CreateExtendableDouble3DArray(hid_t file_id, const char *datasetPath, int stride1, int stride2)
{
    pthread_mutex_lock(&mutexHDF);

    int rank = 3;
    hsize_t initialDimensions[3]  = {stride1, stride2, 0};
    hsize_t maximumDimensions[3]  = {stride1, stride2, H5S_UNLIMITED};

    hid_t dataspace = H5Screate_simple(rank, initialDimensions, maximumDimensions);

    // need chunking for extendable dataset
    hsize_t chunkDimensions[3] = {stride1, stride2, 1};
    hid_t datasetCreationProperty = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(datasetCreationProperty, rank, chunkDimensions);

    hid_t dataset = H5Dcreate2(file_id, datasetPath, H5T_NATIVE_DOUBLE, dataspace,
                               H5P_DEFAULT, datasetCreationProperty, H5P_DEFAULT);

    H5Dclose(dataset);
    H5Sclose(dataspace);

    pthread_mutex_unlock(&mutexHDF);
}

void saveTotalWalltime(const double timeInSeconds) {
    hsize_t dims[1] = {1};

    pthread_mutex_lock(&mutexHDF);
    H5LTmake_dataset(file_id, "/totalWallTimeInSec", 1, dims, H5T_NATIVE_DOUBLE, &timeInSeconds);
    pthread_mutex_unlock(&mutexHDF);
}

void saveLocalOptimizerResults(datapath path, double finalNegLogLikelihood, double masterTime, int exitStatus)
{
    hsize_t dims[1] = {1};

    char fullPath[100];

    pthread_mutex_lock(&mutexHDF);

    sprintf(fullPath, "/crossvalidations/%d/multistarts/%d/finalNegLogLikelihood", path.idxMultiStart, path.idxLocalOptimization);
    H5LTmake_dataset(file_id, fullPath, 1, dims, H5T_NATIVE_DOUBLE, &finalNegLogLikelihood);

    sprintf(fullPath, "/crossvalidations/%d/multistarts/%d/timeOnMaster", path.idxMultiStart, path.idxLocalOptimization);
    H5LTmake_dataset(file_id, fullPath, 1, dims, H5T_NATIVE_DOUBLE, &masterTime);

    sprintf(fullPath, "/crossvalidations/%d/multistarts/%d/exitStatus", path.idxMultiStart, path.idxLocalOptimization);
    H5LTmake_dataset(file_id, fullPath, 1, dims, H5T_NATIVE_INT, &exitStatus);

    pthread_mutex_unlock(&mutexHDF);

    flushResultWriter();
}
