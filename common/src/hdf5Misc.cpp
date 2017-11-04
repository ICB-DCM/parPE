#include "hdf5Misc.h"
#include "logging.h"
#include <assert.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/stat.h>
#include <H5Cpp.h>

namespace parpe {

// mutex for **ALL** HDF5 library calls; read and write; any file(?)
static mutexHdfType mutexHdf;

void initHDF5Mutex() {
    // TODO: check if still required
    H5dont_atexit();
}

void hdf5LockMutex() { mutexHdf.lock(); }

void hdf5UnlockMutex() { mutexHdf.unlock(); }

std::unique_lock<mutexHdfType> hdf5MutexGetLock()
{
    return std::unique_lock<mutexHdfType>(mutexHdf);
}

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
    std::lock_guard<mutexHdfType> lock(mutexHdf);
    bool exists = H5Lexists(file_id, datasetName, H5P_DEFAULT) > 0;

    return exists;
}

bool hdf5GroupExists(hid_t file_id, const char *groupName) {
    std::lock_guard<mutexHdfType> lock(mutexHdf);

    // switch off error handler, check existance and reenable
    H5_SAVE_ERROR_HANDLER;

    herr_t status = H5Gget_objinfo(file_id, groupName, 0, NULL);

    H5_RESTORE_ERROR_HANDLER;

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

    std::lock_guard<mutexHdfType> lock(mutexHdf);

    if (recursively) {
        groupCreationPropertyList = H5Pcreate(H5P_LINK_CREATE);
        H5Pset_create_intermediate_group(groupCreationPropertyList, 1);
    }

    hid_t group = H5Gcreate(file_id, groupPath, groupCreationPropertyList,
                            H5P_DEFAULT, H5P_DEFAULT);
    if (group < 0)
        throw(HDF5Exception("Failed to create group in hdf5CreateGroup: %s", groupPath));
    H5Gclose(group);

}

void hdf5CreateExtendableDouble2DArray(hid_t file_id, const char *datasetPath,
                                       hsize_t stride) {
    int rank = 2;
    hsize_t initialDimensions[2] = {stride, 0};
    hsize_t maximumDimensions[2] = {stride, H5S_UNLIMITED};

    std::lock_guard<mutexHdfType> lock(mutexHdf);

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

}

void hdf5Extend2ndDimensionAndWriteToDouble2DArray(hid_t file_id,
                                                   const char *datasetPath,
                                                   const double *buffer) {
    std::lock_guard<mutexHdfType> lock(mutexHdf);

    hid_t dataset = H5Dopen2(file_id, datasetPath, H5P_DEFAULT);
    if (dataset < 0) {
        throw HDF5Exception("Failed to open dataset %s in hdf5Extend2ndDimensionAndWriteToDouble2DArray", datasetPath);
    }

    // extend
    hid_t filespace = H5Dget_space(dataset);
    int rank = H5Sget_simple_extent_ndims(filespace);
    if (rank != 2) {
        throw HDF5Exception("Failed to write data in hdf5Extend2ndDimensionAndWriteToDouble2DArray: not of rank 2 (%d) when writing %s",
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
            throw HDF5Exception("Failed to write data in "
                  "hdf5Extend2ndDimensionAndWriteToDouble2DArray");

        H5Sclose(memspace);
    }

    H5Sclose(filespace);
    H5Dclose(dataset);
}

void hdf5Extend3rdDimensionAndWriteToDouble3DArray(hid_t file_id,
                                                   const char *datasetPath,
                                                   const double *buffer) {
    std::lock_guard<mutexHdfType> lock(mutexHdf);

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
        throw HDF5Exception("Failed to write data in "
              "hdf5Extend3rdDimensionAndWriteToDouble3DArray");

    H5Dclose(dataset);
    H5Sclose(filespace);
    H5Sclose(memspace);
}

void hdf5CreateOrExtendAndWriteToDouble2DArray(hid_t file_id,
                                               const char *parentPath,
                                               const char *datasetName,
                                               const double *buffer,
                                               hsize_t stride) {
    hdf5EnsureGroupExists(file_id, parentPath);

    std::string fullDatasetPath = std::string(parentPath) + datasetName;

    if (!hdf5DatasetExists(file_id, fullDatasetPath.c_str())) {
        hdf5CreateExtendableDouble2DArray(file_id, fullDatasetPath.c_str(), stride);
    }

    hdf5Extend2ndDimensionAndWriteToDouble2DArray(file_id, fullDatasetPath.c_str(),
                                                  buffer);
}

void hdf5CreateOrExtendAndWriteToDouble3DArray(hid_t file_id,
                                               const char *parentPath,
                                               const char *datasetName,
                                               const double *buffer,
                                               hsize_t stride1, hsize_t stride2) {
    hdf5EnsureGroupExists(file_id, parentPath);

    std::string fullDatasetPath = std::string(parentPath) + datasetName;

    if (!hdf5DatasetExists(file_id, fullDatasetPath.c_str())) {
        hdf5CreateExtendableDouble3DArray(file_id, fullDatasetPath.c_str(), stride1,
                                          stride2);
    }

    hdf5Extend3rdDimensionAndWriteToDouble3DArray(file_id, fullDatasetPath.c_str(),
                                                  buffer);

}

void hdf5CreateOrExtendAndWriteToInt2DArray(hid_t file_id,
                                            const char *parentPath,
                                            const char *datasetName,
                                            const int *buffer, hsize_t stride) {
    std::lock_guard<mutexHdfType> lock(mutexHdf);

    hdf5EnsureGroupExists(file_id, parentPath);

    auto fullDatasetPath = std::string(parentPath) + datasetName;

    if (!hdf5DatasetExists(file_id, fullDatasetPath.c_str())) {
        hdf5CreateExtendableInt2DArray(file_id, fullDatasetPath.c_str(), stride);
    }

    hdf5Extend2ndDimensionAndWriteToInt2DArray(file_id, fullDatasetPath.c_str(),
                                               buffer);
}

void hdf5Extend2ndDimensionAndWriteToInt2DArray(hid_t file_id,
                                                const char *datasetPath,
                                                const int *buffer) {
    std::lock_guard<mutexHdfType> lock(mutexHdf);

    hid_t dataset = H5Dopen2(file_id, datasetPath, H5P_DEFAULT);
    if (dataset < 0)
        throw HDF5Exception("Unable to open dataset %s", datasetPath);

    // extend
    hid_t filespace = H5Dget_space(dataset);
    int rank = H5Sget_simple_extent_ndims(filespace);
    if(rank != 2)
        throw HDF5Exception("Only works for 2D arrays!");

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
        throw HDF5Exception("Error writing data in "
              "hdf5Extend2ndDimensionAndWriteToInt2DArray.");

    H5Sclose(filespace);
    H5Sclose(memspace);

    H5Dclose(dataset);

}

void hdf5CreateExtendableInt2DArray(hid_t file_id, const char *datasetPath,
                                    hsize_t stride) {
    std::lock_guard<mutexHdfType> lock(mutexHdf);

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
}

void hdf5CreateExtendableDouble3DArray(hid_t file_id, const char *datasetPath,
                                       hsize_t stride1, hsize_t stride2) {

    int rank = 3;
    hsize_t initialDimensions[3] = {stride1, stride2, 0};
    hsize_t maximumDimensions[3] = {stride1, stride2, H5S_UNLIMITED};

    std::lock_guard<mutexHdfType> lock(mutexHdf);

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
}

int hdf5Read2DDoubleHyperslab(hid_t file_id, const char *path, hsize_t size0,
                              hsize_t size1, hsize_t offset0, hsize_t offset1,
                              double *buffer) {
    std::lock_guard<mutexHdfType> lock(mutexHdf);

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
    std::lock_guard<mutexHdfType> lock(mutexHdf);

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

int hdf5AttributeExists(hid_t fileId, const char *datasetPath,
                        const char *attributeName) {
    std::lock_guard<mutexHdfType> lock(mutexHdf);

    H5_SAVE_ERROR_HANDLER;
    if (H5Lexists(fileId, datasetPath, H5P_DEFAULT)) {
        hid_t loc = H5Oopen(fileId, datasetPath, H5P_DEFAULT);
        if (loc > 0 && H5LTfind_attribute(loc, attributeName))
            return 1;
    }
    H5_RESTORE_ERROR_HANDLER;
    return 0;
}

void hdf5WriteStringAttribute(hid_t fileId, const char *datasetPath,
                             const char *attributeName,
                             const char *attributeValue) {
    std::lock_guard<mutexHdfType> lock(mutexHdf);

    int ret = H5LTset_attribute_string(fileId, datasetPath, attributeName,
                                    attributeValue);
    if(ret < 0)
        throw HDF5Exception("Unable to write attribute %s on %s",
                            datasetPath, attributeName);
}

hid_t hdf5OpenFile(const char *filename, bool overwrite)
{

    if (!overwrite) {
        struct stat st = {0};
        bool fileExists = stat(filename, &st) == 0;

        if(fileExists)
            throw HDF5Exception("Result file exists");
    }

    std::lock_guard<mutexHdfType> lock(mutexHdf);

    H5_SAVE_ERROR_HANDLER;
    hid_t file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    if (file_id < 0) {
        H5Eprint(H5E_DEFAULT, stderr);
        //throw HDF5Exception();
    }
    H5_RESTORE_ERROR_HANDLER;

    return file_id;
}


void hdf5GetDatasetDimensions(hid_t file_id, const char *path, hsize_t nDimsExpected, int *d1, int *d2, int *d3, int *d4)
{
    std::lock_guard<mutexHdfType> lock(mutexHdf);
    H5_SAVE_ERROR_HANDLER;

    auto dataspace = H5::H5File(file_id).openDataSet(path).getSpace();

    const int nDimsActual = dataspace.getSimpleExtentNdims();
    if(nDimsActual != (signed)nDimsExpected)
        throw(HDF5Exception("Dataset rank does not match nDims argument"));

    hsize_t dims[nDimsExpected];
    dataspace.getSimpleExtentDims(dims, nullptr);

    if(nDimsExpected > 0 && d1)
        *d1 = dims[0];
    if(nDimsExpected > 1 && d2)
        *d2 = dims[1];
    if(nDimsExpected > 2 && d3)
        *d3 = dims[2];
    if(nDimsExpected > 3 && d4)
        *d4 = dims[3];

    H5_RESTORE_ERROR_HANDLER;
}

HDF5Exception::HDF5Exception(const char *format, ...) {
    va_list argptr;
    va_start(argptr,format);
    size_t needed = vsnprintf(nullptr, 0, format, argptr) + 1;
    char buf[needed];

    va_start(argptr,format);
    vsprintf(buf, format, argptr);
    va_end(argptr);

    msg = buf;
}

} // namespace parpe
