#include "hdf5Misc.h"

#include "logging.h"
#include "misc.h"
#include <assert.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>

// mutex for **ALL** HDF5 library calls; read and write; any file(?)
static pthread_mutex_t mutexHDF;

static char *myStringCat(const char *first, const char *second);

void initHDF5Mutex() {
    pthread_mutexattr_t attr;
    pthread_mutexattr_init(&attr);
    pthread_mutexattr_settype(&attr, PTHREAD_MUTEX_RECURSIVE);
    pthread_mutex_init(&mutexHDF, &attr);
    pthread_mutexattr_destroy(&attr);

    H5dont_atexit();
}

void hdf5LockMutex() {
    pthread_mutex_lock(&mutexHDF);
}

void hdf5UnlockMutex() {
    pthread_mutex_unlock(&mutexHDF);
}

void destroyHDF5Mutex() {
    pthread_mutex_destroy(&mutexHDF);
}

herr_t hdf5ErrorStackWalker_cb(unsigned int n, const H5E_error_t *err_desc, void *client_data)
{
    assert (err_desc);
    const int		indent = 2;

    const char *maj_str = H5Eget_major(err_desc->maj_num);
    const char *min_str = H5Eget_minor(err_desc->min_num);

    logmessage(LOGLVL_CRITICAL, "%*s#%03d: %s line %u in %s(): %s",
         indent, "", n, err_desc->file_name, err_desc->line,
         err_desc->func_name, err_desc->desc);
    logmessage(LOGLVL_CRITICAL, "%*smajor(%02d): %s",
         indent*2, "", err_desc->maj_num, maj_str);
    logmessage(LOGLVL_CRITICAL, "%*sminor(%02d): %s",
         indent*2, "", err_desc->min_num, min_str);

    return 0;
}

bool hdf5DatasetExists(hid_t file_id, const char *datasetName)
{
    hdf5LockMutex();
    bool exists = H5Lexists(file_id, datasetName, H5P_DEFAULT) > 0;
    hdf5UnlockMutex();

    return exists;
}

bool hdf5GroupExists(hid_t file_id, const char *groupName)
{
    hdf5LockMutex();

    // switch off error handler, check existance and reenable
    H5_SAVE_ERROR_HANDLER;

    herr_t status = H5Gget_objinfo(file_id, groupName, 0, NULL);

    H5_RESTORE_ERROR_HANDLER;

    hdf5UnlockMutex();

    return status >= 0;
}

void hdf5CreateGroup(hid_t file_id, const char *groupPath, bool recursively)
{
    hid_t groupCreationPropertyList = H5P_DEFAULT;

    hdf5LockMutex();

    if(recursively) {
        groupCreationPropertyList = H5Pcreate(H5P_LINK_CREATE);
        H5Pset_create_intermediate_group (groupCreationPropertyList, 1);
    }

    hid_t group = H5Gcreate(file_id, groupPath, groupCreationPropertyList, H5P_DEFAULT, H5P_DEFAULT);
    if(group < 0)
        error("Failed to create group in hdf5CreateGroup");
    H5Gclose(group);

    hdf5UnlockMutex();
}

void hdf5CreateExtendableDouble2DArray(hid_t file_id, const char *datasetPath, int stride)
{
    int rank = 2;
    hsize_t initialDimensions[2]  = {stride, 0};
    hsize_t maximumDimensions[2]  = {stride, H5S_UNLIMITED};

    hdf5LockMutex();

    hid_t dataspace = H5Screate_simple(rank, initialDimensions, maximumDimensions);

    // need chunking for extendable dataset
    hsize_t chunkDimensions[2] = {stride, 1};
    hid_t datasetCreationProperty = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(datasetCreationProperty, rank, chunkDimensions);

    hid_t dataset = H5Dcreate2(file_id, datasetPath, H5T_NATIVE_DOUBLE, dataspace,
                               H5P_DEFAULT, datasetCreationProperty, H5P_DEFAULT);

    H5Dclose(dataset);
    H5Sclose(dataspace);

    hdf5UnlockMutex();
}

void hdf5Extend2ndDimensionAndWriteToDouble2DArray(hid_t file_id, const char *datasetPath, const double *buffer)
{
    hdf5LockMutex();

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

    hdf5UnlockMutex();
}

void hdf5Extend3rdDimensionAndWriteToDouble3DArray(hid_t file_id, const char *datasetPath, const double *buffer)
{
    hdf5LockMutex();

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

    hdf5UnlockMutex();
}



void hdf5CreateOrExtendAndWriteToDouble2DArray(hid_t file_id, const char *parentPath, const char *datasetName, const double *buffer, int stride)
{
    hdf5LockMutex();

    if(!hdf5GroupExists(file_id, parentPath)) {
        hdf5CreateGroup(file_id, parentPath, true);
    }

    char *fullDatasetPath = myStringCat(parentPath, datasetName);

    if(!hdf5DatasetExists(file_id, fullDatasetPath)) {
        hdf5CreateExtendableDouble2DArray(file_id, fullDatasetPath, stride);
    }

    hdf5Extend2ndDimensionAndWriteToDouble2DArray(file_id, fullDatasetPath, buffer);

    free(fullDatasetPath);

    hdf5UnlockMutex();
}

void hdf5CreateOrExtendAndWriteToDouble3DArray(hid_t file_id, const char *parentPath, const char *datasetName, const double *buffer, int stride1, int stride2)
{
    hdf5LockMutex();

    if(!hdf5GroupExists(file_id, parentPath)) {
        hdf5CreateGroup(file_id, parentPath, true);
    }

    char *fullDatasetPath = myStringCat(parentPath, datasetName);

    if(!hdf5DatasetExists(file_id, fullDatasetPath)) {
        hdf5CreateExtendableDouble3DArray(file_id, fullDatasetPath, stride1, stride2);
    }

    hdf5Extend3rdDimensionAndWriteToDouble3DArray(file_id, fullDatasetPath, buffer);

    hdf5UnlockMutex();

    free(fullDatasetPath);
}


void hdf5CreateOrExtendAndWriteToInt2DArray(hid_t file_id, const char *parentPath, const char *datasetName, const int *buffer, int stride)
{
    hdf5LockMutex();

    if(!hdf5GroupExists(file_id, parentPath)) {
        hdf5CreateGroup(file_id, parentPath, true);
    }

    char *fullDatasetPath = myStringCat(parentPath, datasetName);

    if(!hdf5DatasetExists(file_id, fullDatasetPath)) {
        hdf5CreateExtendableInt2DArray(file_id, fullDatasetPath, stride);
    }

    hdf5Extend2ndDimensionAndWriteToInt2DArray(file_id, fullDatasetPath, buffer);

    hdf5UnlockMutex();

    free(fullDatasetPath);
}

void hdf5Extend2ndDimensionAndWriteToInt2DArray(hid_t file_id, const char *datasetPath, const int *buffer)
{
    hdf5LockMutex();

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

    hdf5UnlockMutex();
}

void hdf5CreateExtendableInt2DArray(hid_t file_id, const char *datasetPath, int stride)
{
    hdf5LockMutex();

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

    hdf5UnlockMutex();
}

void hdf5CreateExtendableDouble3DArray(hid_t file_id, const char *datasetPath, int stride1, int stride2)
{
    hdf5LockMutex();

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

    hdf5UnlockMutex();
}

char *myStringCat(const char *first, const char *second)
{
    char *concatenation = malloc(sizeof(char) * (strlen(first) + strlen(second) + 1));
    strcpy(concatenation, first);
    strcat(concatenation, second);
    return concatenation;
}
