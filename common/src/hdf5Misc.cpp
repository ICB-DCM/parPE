#include "hdf5Misc.h"

#include "logging.h"
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/stat.h>

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

void hdf5LockMutex() { pthread_mutex_lock(&mutexHDF); }

void hdf5UnlockMutex() { pthread_mutex_unlock(&mutexHDF); }

void destroyHDF5Mutex() { pthread_mutex_destroy(&mutexHDF); }

herr_t hdf5ErrorStackWalker_cb(unsigned int n, const H5E_error_t *err_desc,
                               void *client_data) {
    assert(err_desc);
    const int indent = 2;

    char *maj_str = H5Eget_major(err_desc->maj_num);
    char *min_str = H5Eget_minor(err_desc->min_num);

    logmessage(LOGLVL_CRITICAL, "%*s#%03d: %s line %u in %s(): %s", indent, "",
               n, err_desc->file_name, err_desc->line, err_desc->func_name,
               err_desc->desc);
    logmessage(LOGLVL_CRITICAL, "%*smajor(%02d): %s", indent * 2, "",
               err_desc->maj_num, maj_str);
    logmessage(LOGLVL_CRITICAL, "%*sminor(%02d): %s", indent * 2, "",
               err_desc->min_num, min_str);
    free(maj_str);
    free(min_str);
    return 0;
}

bool hdf5DatasetExists(hid_t file_id, const char *datasetName) {
    hdf5LockMutex();
    bool exists = H5Lexists(file_id, datasetName, H5P_DEFAULT) > 0;
    hdf5UnlockMutex();

    return exists;
}

bool hdf5GroupExists(hid_t file_id, const char *groupName) {
    hdf5LockMutex();

    // switch off error handler, check existance and reenable
    H5_SAVE_ERROR_HANDLER;

    herr_t status = H5Gget_objinfo(file_id, groupName, 0, NULL);

    H5_RESTORE_ERROR_HANDLER;

    hdf5UnlockMutex();

    return status >= 0;
}

bool hdf5EnsureGroupExists(hid_t file_id, const char *groupName) {
    if (!hdf5GroupExists(file_id, groupName)) {
        hdf5CreateGroup(file_id, groupName, true);
    }

    return 0;
}

void hdf5CreateGroup(hid_t file_id, const char *groupPath, bool recursively) {
    hid_t groupCreationPropertyList = H5P_DEFAULT;

    hdf5LockMutex();

    if (recursively) {
        groupCreationPropertyList = H5Pcreate(H5P_LINK_CREATE);
        H5Pset_create_intermediate_group(groupCreationPropertyList, 1);
    }

    hid_t group = H5Gcreate(file_id, groupPath, groupCreationPropertyList,
                            H5P_DEFAULT, H5P_DEFAULT);
    if (group < 0)
        logmessage(LOGLVL_ERROR,
                   "Failed to create group in hdf5CreateGroup: %s", groupPath);
    H5Gclose(group);

    hdf5UnlockMutex();
}

void hdf5CreateExtendableDouble2DArray(hid_t file_id, const char *datasetPath,
                                       hsize_t stride) {
    int rank = 2;
    hsize_t initialDimensions[2] = {stride, 0};
    hsize_t maximumDimensions[2] = {stride, H5S_UNLIMITED};

    hdf5LockMutex();

    hid_t dataspace =
        H5Screate_simple(rank, initialDimensions, maximumDimensions);

    // need chunking for extendable dataset
    hsize_t chunkDimensions[2] = {stride, 1};
    hid_t datasetCreationProperty = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(datasetCreationProperty, rank, chunkDimensions);

    hid_t dataset =
        H5Dcreate2(file_id, datasetPath, H5T_NATIVE_DOUBLE, dataspace,
                   H5P_DEFAULT, datasetCreationProperty, H5P_DEFAULT);

    H5Dclose(dataset);
    H5Sclose(dataspace);

    hdf5UnlockMutex();
}

void hdf5Extend2ndDimensionAndWriteToDouble2DArray(hid_t file_id,
                                                   const char *datasetPath,
                                                   const double *buffer) {
    hdf5LockMutex();

    hid_t dataset = H5Dopen2(file_id, datasetPath, H5P_DEFAULT);
    if (dataset < 0) {
        logmessage(LOGLVL_CRITICAL, "Failed to open dataset %s in "
                                    "hdf5Extend2ndDimensionAndWriteToDouble2DAr"
                                    "ray",
                   datasetPath);
        hdf5UnlockMutex();
        return;
    }

    // extend
    hid_t filespace = H5Dget_space(dataset);
    int rank = H5Sget_simple_extent_ndims(filespace);
    if (rank != 2) {
        logmessage(LOGLVL_CRITICAL, "Failed to write data in "
                                    "hdf5Extend2ndDimensionAndWriteToDouble2DAr"
                                    "ray: not of rank 2 (%d) when writing %s",
                   rank, datasetPath);
    } else {

        hsize_t currentDimensions[2];
        H5Sget_simple_extent_dims(filespace, currentDimensions, NULL);

        hsize_t newDimensions[2] = {currentDimensions[0], currentDimensions[1] + 1};
        herr_t status = H5Dset_extent(dataset, newDimensions);

        filespace = H5Dget_space(dataset);
        hsize_t offset[2] = {0, currentDimensions[1]};
        hsize_t slabsize[2] = {currentDimensions[0], 1};

        status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL,
                                     slabsize, NULL);

        hid_t memspace = H5Screate_simple(rank, slabsize, NULL);

        status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, memspace, filespace,
                          H5P_DEFAULT, buffer);

        if (status < 0)
            error("Failed to write data in "
                  "hdf5Extend2ndDimensionAndWriteToDouble2DArray");

        H5Sclose(memspace);
    }

    H5Sclose(filespace);
    H5Dclose(dataset);
    hdf5UnlockMutex();
}

void hdf5Extend3rdDimensionAndWriteToDouble3DArray(hid_t file_id,
                                                   const char *datasetPath,
                                                   const double *buffer) {
    hdf5LockMutex();

    hid_t dataset = H5Dopen2(file_id, datasetPath, H5P_DEFAULT);

    // extend
    hid_t filespace = H5Dget_space(dataset);
    int rank = H5Sget_simple_extent_ndims(filespace);
    assert(rank == 3 && "Only works for 3D arrays!");

    hsize_t currentDimensions[3];
    H5Sget_simple_extent_dims(filespace, currentDimensions, NULL);

    hsize_t newDimensions[3] = {currentDimensions[0], currentDimensions[1],
                                currentDimensions[2] + 1};
    herr_t status = H5Dset_extent(dataset, newDimensions);

    filespace = H5Dget_space(dataset);
    hsize_t offset[3] = {0, 0, currentDimensions[2]};
    hsize_t slabsize[3] = {currentDimensions[0], currentDimensions[1], 1};

    status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL,
                                 slabsize, NULL);

    hid_t memspace = H5Screate_simple(rank, slabsize, NULL);

    status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, memspace, filespace,
                      H5P_DEFAULT, buffer);

    if (status < 0)
        error("Failed to write data in "
              "hdf5Extend3rdDimensionAndWriteToDouble3DArray");

    H5Dclose(dataset);
    H5Sclose(filespace);
    H5Sclose(memspace);

    hdf5UnlockMutex();
}

void hdf5CreateOrExtendAndWriteToDouble2DArray(hid_t file_id,
                                               const char *parentPath,
                                               const char *datasetName,
                                               const double *buffer,
                                               hsize_t stride) {
    hdf5LockMutex();

    hdf5EnsureGroupExists(file_id, parentPath);

    char *fullDatasetPath = myStringCat(parentPath, datasetName);

    if (!hdf5DatasetExists(file_id, fullDatasetPath)) {
        hdf5CreateExtendableDouble2DArray(file_id, fullDatasetPath, stride);
    }

    hdf5Extend2ndDimensionAndWriteToDouble2DArray(file_id, fullDatasetPath,
                                                  buffer);

    delete[] fullDatasetPath;

    hdf5UnlockMutex();
}

void hdf5CreateOrExtendAndWriteToDouble3DArray(hid_t file_id,
                                               const char *parentPath,
                                               const char *datasetName,
                                               const double *buffer,
                                               hsize_t stride1, hsize_t stride2) {
    hdf5LockMutex();

    hdf5EnsureGroupExists(file_id, parentPath);

    char *fullDatasetPath = myStringCat(parentPath, datasetName);

    if (!hdf5DatasetExists(file_id, fullDatasetPath)) {
        hdf5CreateExtendableDouble3DArray(file_id, fullDatasetPath, stride1,
                                          stride2);
    }

    hdf5Extend3rdDimensionAndWriteToDouble3DArray(file_id, fullDatasetPath,
                                                  buffer);

    hdf5UnlockMutex();

    delete[] fullDatasetPath;
}

void hdf5CreateOrExtendAndWriteToInt2DArray(hid_t file_id,
                                            const char *parentPath,
                                            const char *datasetName,
                                            const int *buffer, hsize_t stride) {
    hdf5LockMutex();

    hdf5EnsureGroupExists(file_id, parentPath);

    char *fullDatasetPath = myStringCat(parentPath, datasetName);

    if (!hdf5DatasetExists(file_id, fullDatasetPath)) {
        hdf5CreateExtendableInt2DArray(file_id, fullDatasetPath, stride);
    }

    hdf5Extend2ndDimensionAndWriteToInt2DArray(file_id, fullDatasetPath,
                                               buffer);

    hdf5UnlockMutex();

    delete[] fullDatasetPath;
}

void hdf5Extend2ndDimensionAndWriteToInt2DArray(hid_t file_id,
                                                const char *datasetPath,
                                                const int *buffer) {
    hdf5LockMutex();

    hid_t dataset = H5Dopen2(file_id, datasetPath, H5P_DEFAULT);
    if (dataset < 0) {
        logmessage(LOGLVL_CRITICAL, "Unable to open dataset %s", datasetPath);
    } else {

        // extend
        hid_t filespace = H5Dget_space(dataset);
        int rank = H5Sget_simple_extent_ndims(filespace);
        assert(rank == 2 && "Only works for 2D arrays!");

        hsize_t currentDimensions[2];
        H5Sget_simple_extent_dims(filespace, currentDimensions, NULL);

        hsize_t newDimensions[2] = {currentDimensions[0], currentDimensions[1] + 1};
        herr_t status = H5Dset_extent(dataset, newDimensions);

        filespace = H5Dget_space(dataset);
        hsize_t offset[2] = {0, currentDimensions[1]};
        hsize_t slabsize[2] = {currentDimensions[0], 1};

        status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL,
                                     slabsize, NULL);

        hid_t memspace = H5Screate_simple(rank, slabsize, NULL);

        status = H5Dwrite(dataset, H5T_NATIVE_INT, memspace, filespace, H5P_DEFAULT,
                          buffer);

        if (status < 0)
            error("Error writing data in "
                  "hdf5Extend2ndDimensionAndWriteToInt2DArray.");

        H5Sclose(filespace);
        H5Sclose(memspace);

    }
    H5Dclose(dataset);

    hdf5UnlockMutex();
}

void hdf5CreateExtendableInt2DArray(hid_t file_id, const char *datasetPath,
                                    hsize_t stride) {
    hdf5LockMutex();

    int rank = 2;
    hsize_t initialDimensions[2] = {stride, 0};
    hsize_t maximumDimensions[2] = {stride, H5S_UNLIMITED};

    hid_t dataspace =
        H5Screate_simple(rank, initialDimensions, maximumDimensions);

    // need chunking for extendable dataset
    hsize_t chunkDimensions[2] = {stride, 1};
    hid_t datasetCreationProperty = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(datasetCreationProperty, rank, chunkDimensions);

    assert(H5Tget_size(H5T_NATIVE_INT) == sizeof(int));
    hid_t dataset =
        H5Dcreate2(file_id, datasetPath, H5T_NATIVE_INT, dataspace, H5P_DEFAULT,
                   datasetCreationProperty, H5P_DEFAULT);
    assert(dataset >= 0 && "Unable to open dataset!");

    H5Dclose(dataset);
    H5Sclose(dataspace);

    hdf5UnlockMutex();
}

void hdf5CreateExtendableDouble3DArray(hid_t file_id, const char *datasetPath,
                                       hsize_t stride1, hsize_t stride2) {
    hdf5LockMutex();

    int rank = 3;
    hsize_t initialDimensions[3] = {stride1, stride2, 0};
    hsize_t maximumDimensions[3] = {stride1, stride2, H5S_UNLIMITED};

    hid_t dataspace =
        H5Screate_simple(rank, initialDimensions, maximumDimensions);

    // need chunking for extendable dataset
    hsize_t chunkDimensions[3] = {stride1, stride2, 1};
    hid_t datasetCreationProperty = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(datasetCreationProperty, rank, chunkDimensions);

    hid_t dataset =
        H5Dcreate2(file_id, datasetPath, H5T_NATIVE_DOUBLE, dataspace,
                   H5P_DEFAULT, datasetCreationProperty, H5P_DEFAULT);

    H5Dclose(dataset);
    H5Sclose(dataspace);

    hdf5UnlockMutex();
}

char *myStringCat(const char *first, const char *second) {
    char *concatenation = new char[strlen(first) + strlen(second) + 1];
    strcpy(concatenation, first);
    strcat(concatenation, second);
    return concatenation;
}

int hdf5Read2DDoubleHyperslab(hid_t file_id, const char *path, hsize_t size0,
                              hsize_t size1, hsize_t offset0, hsize_t offset1,
                              double *buffer) {
    hid_t dataset = H5Dopen2(file_id, path, H5P_DEFAULT);
    hid_t dataspace = H5Dget_space(dataset);
    hsize_t offset[] = {offset0, offset1};
    hsize_t count[] = {size0, size1};

    const int ndims = H5Sget_simple_extent_ndims(dataspace);
    assert(ndims == 2 && "Only works for 2D arrays!");
    hsize_t dims[ndims];
    H5Sget_simple_extent_dims(dataspace, dims, NULL);
    // printf("%lld %lld, %lld %lld, %lld %lld\n", dims[0], dims[1], offset0,
    // offset1, size0, size1);
    assert(dims[0] >= offset0 && dims[0] >= size0 &&
           "Offset larger than dataspace dimensions!");
    assert(dims[1] >= offset1 && dims[1] >= size1 &&
           "Offset larger than dataspace dimensions!");

    H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset, NULL, count, NULL);

    hid_t memspace = H5Screate_simple(2, count, NULL);

    H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace, H5P_DEFAULT,
            buffer);

    H5Sclose(dataspace);
    H5Dclose(dataset);

    return 0;
}

int hdf5Read3DDoubleHyperslab(hid_t file_id, const char *path, hsize_t size0,
                              hsize_t size1, hsize_t size2, hsize_t offset0,
                              hsize_t offset1, hsize_t offset2,
                              double *buffer) {
    const int rank = 3;
    hid_t dataset = H5Dopen2(file_id, path, H5P_DEFAULT);
    hid_t dataspace = H5Dget_space(dataset);
    hsize_t offset[] = {offset0, offset1, offset2};
    hsize_t count[] = {size0, size1, size2};

    const int ndims = H5Sget_simple_extent_ndims(dataspace);
    assert(ndims == rank && "Only works for 3D arrays!");
    hsize_t dims[ndims];
    H5Sget_simple_extent_dims(dataspace, dims, NULL);
    assert(dims[0] >= offset0 && dims[0] >= size0 &&
           "Offset larger than dataspace dimensions!");
    assert(dims[1] >= offset1 && dims[1] >= size1 &&
           "Offset larger than dataspace dimensions!");
    assert(dims[2] >= offset2 && dims[2] >= size2 &&
           "Offset larger than dataspace dimensions!");

    H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset, NULL, count, NULL);

    hid_t memspace = H5Screate_simple(rank, count, NULL);

    H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace, H5P_DEFAULT,
            buffer);

    H5Sclose(dataspace);
    H5Dclose(dataset);

    return 0;
}

void hdf5GetDatasetDimensions2D(hid_t file_id, const char *path, int *d1,
                                int *d2) {
    hdf5LockMutex();

    hid_t dataset = H5Dopen2(file_id, path, H5P_DEFAULT);
    hid_t dataspace = H5Dget_space(dataset);

    const int ndims = H5Sget_simple_extent_ndims(dataspace);
    assert(ndims == 2 && "Only works for 2D arrays!");
    hsize_t dims[ndims];
    H5Sget_simple_extent_dims(dataspace, dims, NULL);

    *d1 = dims[0];
    *d2 = dims[1];

    H5Sclose(dataspace);
    H5Dclose(dataset);

    hdf5UnlockMutex();
}

int hdf5AttributeExists(hid_t fileId, const char *datasetPath,
                        const char *attributeName) {
    H5_SAVE_ERROR_HANDLER;
    if (H5Lexists(fileId, datasetPath, H5P_DEFAULT)) {
        hid_t dataset = H5Dopen2(fileId, datasetPath, H5P_DEFAULT);
        if (H5LTfind_attribute(dataset, attributeName))
            return 1;
    }
    H5_RESTORE_ERROR_HANDLER;
    return 0;
}

int hdf5WriteStringAttribute(hid_t fileId, const char *datasetPath,
                             const char *attributeName,
                             const char *attributeValue) {
    return H5LTset_attribute_string(fileId, datasetPath, attributeName,
                                    attributeValue);
}

hid_t hdf5OpenFile(const char *filename, bool overwrite)
{
    if (!overwrite) {
        struct stat st = {0};
        bool fileExists = stat(filename, &st) == 0;

        // throw HDF5Exception("Result file exists");
        return -1;
        assert(!fileExists);
    }

    H5_SAVE_ERROR_HANDLER;
    hid_t file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    if (file_id < 0) {
        H5Eprint(H5E_DEFAULT, stderr);
        //throw HDF5Exception();
    }
    H5_RESTORE_ERROR_HANDLER;

    return file_id;
}

void hdf5GetDatasetDimensions3D(hid_t file_id, const char *path, int *d1, int *d2, int *d3)
{
    hdf5LockMutex();

    hid_t dataset = H5Dopen2(file_id, path, H5P_DEFAULT);
    hid_t dataspace = H5Dget_space(dataset);

    const int ndims = H5Sget_simple_extent_ndims(dataspace);
    assert(ndims == 3 && "Only works for 3D arrays!");
    hsize_t dims[ndims];
    H5Sget_simple_extent_dims(dataspace, dims, NULL);

    *d1 = dims[0];
    *d2 = dims[1];
    *d3 = dims[2];

    H5Sclose(dataspace);
    H5Dclose(dataset);

    hdf5UnlockMutex();
}
