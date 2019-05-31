#include <parpecommon/hdf5Misc.h>
#include <parpecommon/logging.h>
#include <parpecommon/misc.h>

#include <cassert>
#include <cstdlib>
#include <utility>
#include <unistd.h>
#include <sys/stat.h>

#include <H5Tpublic.h>

namespace parpe {

/** mutex for **ALL** HDF5 library calls; read and write; any file(?) */
static mutexHdfType mutexHdf;

void initHDF5Mutex() {
    // TODO: check if still required
    H5dont_atexit();
}

std::unique_lock<mutexHdfType> hdf5MutexGetLock()
{
    return std::unique_lock<mutexHdfType>(mutexHdf);
}

herr_t hdf5ErrorStackWalker_cb(unsigned int n, const H5E_error_t *err_desc,
                               void* /*client_data*/) {
    assert(err_desc);
    const int indent = 2;

    std::unique_ptr<char, decltype(std::free) *>
            maj_str { H5Eget_major(err_desc->maj_num), std::free };
    std::unique_ptr<char, decltype(std::free) *>
            min_str { H5Eget_minor(err_desc->min_num), std::free };

    logmessage(LOGLVL_CRITICAL, "%*s#%03d: %s line %u in %s(): %s", indent, "",
               n, err_desc->file_name, err_desc->line, err_desc->func_name,
               err_desc->desc);
    logmessage(LOGLVL_CRITICAL, "%*smajor(%02d): %s", indent * 2, "",
               err_desc->maj_num, maj_str.get());
    logmessage(LOGLVL_CRITICAL, "%*sminor(%02d): %s", indent * 2, "",
               err_desc->min_num, min_str.get());

    return 0;
}

bool hdf5DatasetExists(hid_t file_id, const char *datasetName)
{
    std::lock_guard<mutexHdfType> lock(mutexHdf);
    bool exists = H5Lexists(file_id, datasetName, H5P_DEFAULT) > 0;

    return exists;
}

bool hdf5GroupExists(hid_t file_id, const char *groupName)
{
    std::lock_guard<mutexHdfType> lock(mutexHdf);

    // switch off error handler, check existance and reenable
    H5_SAVE_ERROR_HANDLER;

    herr_t status = H5Gget_objinfo(file_id, groupName, false, nullptr);

    H5_RESTORE_ERROR_HANDLER;

    return status >= 0;
}

void hdf5EnsureGroupExists(hid_t file_id, const char *groupName) {
    if (!hdf5GroupExists(file_id, groupName)) {
        hdf5CreateGroup(file_id, groupName, true);
    }
}

void hdf5CreateGroup(hid_t file_id, const char *groupPath, bool recursively)
{
    auto groupCreationPropertyList = H5P_DEFAULT;

    std::lock_guard<mutexHdfType> lock(mutexHdf);

    if (recursively) {
        groupCreationPropertyList = H5Pcreate(H5P_LINK_CREATE);
        H5Pset_create_intermediate_group(groupCreationPropertyList, 1);
    }

    auto group = H5Gcreate(file_id, groupPath, groupCreationPropertyList,
                           H5P_DEFAULT, H5P_DEFAULT);
    if (group < 0)
        throw(HDF5Exception("Failed to create group in hdf5CreateGroup: %s",
                            groupPath));

    H5Gclose(group);
}

void hdf5CreateExtendableDouble2DArray(hid_t file_id,
                                       const char *datasetPath,
                                       hsize_t stride)
{
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

    if(dataset < 0)
        throw HDF5Exception("hdf5CreateExtendableDouble2DArray");

    H5Dclose(dataset);
    H5Sclose(dataspace);

}

void hdf5Extend2ndDimensionAndWriteToDouble2DArray(
        hid_t file_id, const char *datasetPath, gsl::span<const double> buffer)
{
    std::lock_guard<mutexHdfType> lock(mutexHdf);

    hid_t dataset = H5Dopen2(file_id, datasetPath, H5P_DEFAULT);
    if (dataset < 0) {
        throw HDF5Exception("Failed to open dataset %s in "
                            "hdf5Extend2ndDimensionAndWriteToDouble2DArray",
                            datasetPath);
    }

    // check rank
    hid_t filespace = H5Dget_space(dataset);
    int rank = H5Sget_simple_extent_ndims(filespace);
    if (rank != 2) {
        H5Sclose(filespace);
        H5Dclose(dataset);
        throw HDF5Exception("Failed to write data in "
                            "hdf5Extend2ndDimensionAndWriteToDouble2DArray: "
                            "not of rank 2 (%d) when writing %s",
                            rank, datasetPath);
    }

    // extend
    hsize_t currentDimensions[2];
    H5Sget_simple_extent_dims(filespace, currentDimensions, nullptr);
    H5Sclose(filespace);

    RELEASE_ASSERT(buffer.size() == currentDimensions[0], "");

    hsize_t newDimensions[2] = {currentDimensions[0], currentDimensions[1] + 1};
    herr_t status = H5Dset_extent(dataset, newDimensions);

    filespace = H5Dget_space(dataset);
    hsize_t offset[2] = {0, currentDimensions[1]};
    hsize_t slabsize[2] = {currentDimensions[0], 1};

    status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, nullptr,
                                 slabsize, nullptr);

    hid_t memspace = H5Screate_simple(rank, slabsize, nullptr);

    status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, memspace, filespace,
                      H5P_DEFAULT, buffer.data());

    H5Sclose(memspace);
    H5Sclose(filespace);
    H5Dclose(dataset);

    if (status < 0)
        throw HDF5Exception("Failed to write data in "
                            "hdf5Extend2ndDimensionAndWriteToDouble2DArray");
}

void hdf5Extend3rdDimensionAndWriteToDouble3DArray(hid_t file_id,
                                                   const char *datasetPath,
                                                   const double *buffer)
{
    std::lock_guard<mutexHdfType> lock(mutexHdf);

    hid_t dataset = H5Dopen2(file_id, datasetPath, H5P_DEFAULT);

    // extend
    hid_t filespace = H5Dget_space(dataset);
    int rank = H5Sget_simple_extent_ndims(filespace);
    assert(rank == 3 && "Only works for 3D arrays!");

    hsize_t currentDimensions[3];
    H5Sget_simple_extent_dims(filespace, currentDimensions, nullptr);

    hsize_t newDimensions[3] = {currentDimensions[0], currentDimensions[1],
                                currentDimensions[2] + 1};
    herr_t status = H5Dset_extent(dataset, newDimensions);

    filespace = H5Dget_space(dataset);
    hsize_t offset[3] = {0, 0, currentDimensions[2]};
    hsize_t slabsize[3] = {currentDimensions[0], currentDimensions[1], 1};

    status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, nullptr,
                                 slabsize, nullptr);

    hid_t memspace = H5Screate_simple(rank, slabsize, nullptr);

    status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, memspace, filespace,
                      H5P_DEFAULT, buffer);

    H5Dclose(dataset);
    H5Sclose(filespace);
    H5Sclose(memspace);

    if (status < 0)
        throw HDF5Exception("Failed to write data in "
                            "hdf5Extend3rdDimensionAndWriteToDouble3DArray");
}

void hdf5CreateOrExtendAndWriteToDouble2DArray(hid_t file_id,
                                               const char *parentPath,
                                               const char *datasetName,
                                               gsl::span<const double> buffer)
{
    hdf5EnsureGroupExists(file_id, parentPath);

    std::string fullDatasetPath = std::string(parentPath) + "/" + datasetName;

    if (!hdf5DatasetExists(file_id, fullDatasetPath)) {
        hdf5CreateExtendableDouble2DArray(
                    file_id, fullDatasetPath.c_str(), buffer.size());
    }

    hdf5Extend2ndDimensionAndWriteToDouble2DArray(
                file_id, fullDatasetPath.c_str(), buffer);
}

void hdf5CreateOrExtendAndWriteToDouble3DArray(hid_t file_id,
                                               const char *parentPath,
                                               const char *datasetName,
                                               gsl::span<const double> buffer,
                                               hsize_t stride1,
                                               hsize_t stride2) {
    hdf5EnsureGroupExists(file_id, parentPath);

    std::string fullDatasetPath = std::string(parentPath) + "/" + datasetName;

    if (!hdf5DatasetExists(file_id, fullDatasetPath)) {
        hdf5CreateExtendableDouble3DArray(
                    file_id, fullDatasetPath.c_str(), stride1, stride2);
    }

    hdf5Extend3rdDimensionAndWriteToDouble3DArray(file_id,
                                                  fullDatasetPath.c_str(),
                                                  buffer.data());

}

void hdf5CreateOrExtendAndWriteToInt2DArray(hid_t file_id,
                                            const char *parentPath,
                                            const char *datasetName,
                                            gsl::span<const int> buffer)
{
    std::lock_guard<mutexHdfType> lock(mutexHdf);

    hdf5EnsureGroupExists(file_id, parentPath);

    auto fullDatasetPath = std::string(parentPath) + "/" + datasetName;

    if (!hdf5DatasetExists(file_id, fullDatasetPath)) {
        hdf5CreateExtendableInt2DArray(
                    file_id, fullDatasetPath.c_str(), buffer.size());
    }

    hdf5Extend2ndDimensionAndWriteToInt2DArray(file_id, fullDatasetPath.c_str(),
                                               buffer);
}

void hdf5Extend2ndDimensionAndWriteToInt2DArray(hid_t file_id,
                                                const char *datasetPath,
                                                gsl::span<const int> buffer)
{
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
    H5Sget_simple_extent_dims(filespace, currentDimensions, nullptr);
    RELEASE_ASSERT(buffer.size() == currentDimensions[0], "");

    hsize_t newDimensions[2] = {currentDimensions[0], currentDimensions[1] + 1};
    herr_t status = H5Dset_extent(dataset, newDimensions);

    filespace = H5Dget_space(dataset);
    hsize_t offset[2] = {0, currentDimensions[1]};
    hsize_t slabsize[2] = {currentDimensions[0], 1};

    status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, nullptr,
                                 slabsize, nullptr);

    hid_t memspace = H5Screate_simple(rank, slabsize, nullptr);

    status = H5Dwrite(dataset, H5T_NATIVE_INT, memspace, filespace, H5P_DEFAULT,
                      buffer.data());
    H5Sclose(filespace);
    H5Sclose(memspace);

    H5Dclose(dataset);

    if (status < 0)
        throw HDF5Exception("Error writing data in "
                            "hdf5Extend2ndDimensionAndWriteToInt2DArray.");
}

void hdf5CreateExtendableInt2DArray(hid_t file_id,
                                    const char *datasetPath,
                                    hsize_t stride)
{
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
            H5Dcreate2(file_id, datasetPath, H5T_NATIVE_INT, dataspace,
                       H5P_DEFAULT, datasetCreationProperty, H5P_DEFAULT);

    if(dataset < 0)
        throw HDF5Exception("hdf5CreateExtendableInt2DArray");

    H5Dclose(dataset);
    H5Sclose(dataspace);
}

void hdf5CreateExtendableDouble3DArray(hid_t file_id,
                                       const char *datasetPath,
                                       hsize_t stride1,
                                       hsize_t stride2)
{

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

    if(dataset < 0)
        throw HDF5Exception("hdf5CreateExtendableDouble3DArray");

    H5Dclose(dataset);
    H5Sclose(dataspace);
}

int hdf5Read2DDoubleHyperslab(hid_t file_id,
                              const char *path,
                              hsize_t size0,
                              hsize_t size1,
                              hsize_t offset0,
                              hsize_t offset1,
                              gsl::span<double> buffer)
{
    RELEASE_ASSERT(buffer.size() == size0 * size1, "");

    std::lock_guard<mutexHdfType> lock(mutexHdf);

    hid_t dataset = H5Dopen2(file_id, path, H5P_DEFAULT);
    hid_t dataspace = H5Dget_space(dataset);
    hsize_t offset[] = {offset0, offset1};
    hsize_t count[] = {size0, size1};

    const int ndims = H5Sget_simple_extent_ndims(dataspace);
    RELEASE_ASSERT(ndims == 2 && "Only works for 2D arrays!", "");
    hsize_t dims[ndims];
    H5Sget_simple_extent_dims(dataspace, dims, nullptr);
    // printf("%lld %lld, %lld %lld, %lld %lld\n", dims[0], dims[1], offset0,
    // offset1, size0, size1);
    assert(dims[0] >= offset0 && dims[0] >= size0 &&
            "Offset larger than dataspace dimensions!");
    assert(dims[1] >= offset1 && dims[1] >= size1 &&
            "Offset larger than dataspace dimensions!");

    H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset, nullptr,
                        count, nullptr);

    hid_t memspace = H5Screate_simple(2, count, nullptr);

    H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace, H5P_DEFAULT,
            buffer.data());

    H5Sclose(dataspace);
    H5Dclose(dataset);

    return 0;
}

std::vector<int> hdf5Read1DIntegerHyperslab(H5::H5File const& file,
                                            std::string const& path,
                                            hsize_t count,
                                            hsize_t offset)
{
    std::lock_guard<mutexHdfType> lock(mutexHdf);

    H5::DataSet dataset = file.openDataSet(path);
    H5::DataSpace filespace = dataset.getSpace();

    const int ndims = filespace.getSimpleExtentNdims();
    assert(ndims == 1 && "Only works for 1D arrays!");
    hsize_t length;
    filespace.getSimpleExtentDims(&length);

    assert(length >= offset && "Offset larger than dataspace dimensions!");

    filespace.selectHyperslab(H5S_SELECT_SET, &count, &offset);

    H5::DataSpace memspace(1, &count);
    std::vector<int> buffer(count);

    dataset.read(buffer.data(), H5::PredType::NATIVE_INT, memspace, filespace);

    return buffer;
}

std::vector<int> hdf5Read2DIntegerHyperslab(const H5::H5File &file,
                                            std::string const& path,
                                            hsize_t size0,
                                            hsize_t size1,
                                            hsize_t offset0,
                                            hsize_t offset1)
{
    std::lock_guard<mutexHdfType> lock(mutexHdf);

    H5::DataSet dataset = file.openDataSet(path);
    H5::DataSpace filespace = dataset.getSpace();

    const int ndims = filespace.getSimpleExtentNdims();
    assert(ndims == 2 && "Only works for 2D arrays!");
    hsize_t dims[ndims];
    filespace.getSimpleExtentDims(dims);

    hsize_t offset[] = {offset0, offset1};
    hsize_t count[] = {size0, size1};
    // printf("%lld %lld, %lld %lld, %lld %lld\n", dims[0], dims[1], offset0,
    // offset1, size0, size1);
    assert(dims[0] >= offset0 && dims[0] >= size0 &&
            "Offset larger than dataspace dimensions!");
    assert(dims[1] >= offset1 && dims[1] >= size1 &&
            "Offset larger than dataspace dimensions!");

    filespace.selectHyperslab(H5S_SELECT_SET, count, offset);

    H5::DataSpace memspace(2, count);
    std::vector<int> buffer(size0 * size1);

    dataset.read(buffer.data(), H5::PredType::NATIVE_INT, memspace, filespace);

    return buffer;
}

int hdf5Read3DDoubleHyperslab(hid_t file_id,
                              const char *path,
                              hsize_t size0,
                              hsize_t size1,
                              hsize_t size2,
                              hsize_t offset0,
                              hsize_t offset1,
                              hsize_t offset2,
                              double *buffer)
{
    std::lock_guard<mutexHdfType> lock(mutexHdf);

    const int rank = 3;
    hid_t dataset = H5Dopen2(file_id, path, H5P_DEFAULT);
    hid_t dataspace = H5Dget_space(dataset);
    hsize_t offset[] = {offset0, offset1, offset2};
    hsize_t count[] = {size0, size1, size2};

    const int ndims = H5Sget_simple_extent_ndims(dataspace);
    assert(ndims == rank && "Only works for 3D arrays!");
    hsize_t dims[ndims];
    H5Sget_simple_extent_dims(dataspace, dims, nullptr);
    assert(dims[0] >= offset0 && dims[0] >= size0 &&
            "Offset larger than dataspace dimensions!");
    assert(dims[1] >= offset1 && dims[1] >= size1 &&
            "Offset larger than dataspace dimensions!");
    assert(dims[2] >= offset2 && dims[2] >= size2 &&
            "Offset larger than dataspace dimensions!");

    H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset, nullptr,
                        count, nullptr);

    hid_t memspace = H5Screate_simple(rank, count, nullptr);

    H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace, H5P_DEFAULT,
            buffer);

    H5Sclose(dataspace);
    H5Dclose(dataset);

    return 0;
}


std::vector<double> hdf5Get3DDoubleHyperslab(hid_t file_id, const char *path,
                                             hsize_t size0, hsize_t size1,
                                             hsize_t size2,
                                             hsize_t offset0, hsize_t offset1,
                                             hsize_t offset2)
{
    std::vector<double> buffer(size0 * size1 * size2);
    hdf5Read3DDoubleHyperslab(file_id, path, size0, size1, size2,
                              offset0, offset1, offset2, buffer.data());
    return buffer;
}


bool hdf5AttributeExists(hid_t fileId,
                         const char *datasetPath,
                         const char *attributeName)
{
    std::lock_guard<mutexHdfType> lock(mutexHdf);

    int exists = false;

    H5_SAVE_ERROR_HANDLER;

    hid_t loc = H5Oopen(fileId, datasetPath, H5P_DEFAULT);
    if (loc >= 0) {
        exists = H5LTfind_attribute(loc, attributeName);
        H5Oclose(loc);
    }

    H5_RESTORE_ERROR_HANDLER;

    return exists;
}

void hdf5WriteStringAttribute(hid_t fileId,
                              const char *datasetPath,
                              const char *attributeName,
                              const char *attributeValue)
{
    std::lock_guard<mutexHdfType> lock(mutexHdf);

    int ret = H5LTset_attribute_string(fileId, datasetPath, attributeName,
                                       attributeValue);
    if(ret < 0)
        throw HDF5Exception("Unable to write attribute %s on %s",
                            datasetPath, attributeName);
}

hid_t hdf5CreateFile(const char *filename,
                     bool overwrite)
{
    // Create parent folders
    mkpathConstChar(filename, 0755);

    std::lock_guard<mutexHdfType> lock(mutexHdf);

    if (!overwrite) {
        struct stat st = {};
        bool fileExists = stat(filename, &st) == 0;

        if(fileExists)
            throw HDF5Exception("Result file exists %s", filename);
    }

    H5_SAVE_ERROR_HANDLER;
    hid_t file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    if (file_id < 0) {
        H5Eprint(H5E_DEFAULT, stderr);
        printBacktrace();
        throw HDF5Exception("hdf5CreateFile: Failed to create file %s. "
                            "Is this file opened by another process?",
                            filename);
    }
    H5_RESTORE_ERROR_HANDLER;

    return file_id;
}


void hdf5GetDatasetDimensions(hid_t file_id,
                              const char *path,
                              hsize_t nDimsExpected,
                              int *d1,
                              int *d2,
                              int *d3,
                              int *d4)
{
    assert(file_id >= 0);

    std::lock_guard<mutexHdfType> lock(mutexHdf);
    H5_SAVE_ERROR_HANDLER;

    auto file = H5::H5File(file_id);
    auto dataset = file.openDataSet(path);
    auto dataspace = dataset.getSpace();

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

HDF5Exception::HDF5Exception(std::string msg) : msg(std::move(msg)) {
    stackTrace = getBacktrace(20);
}

HDF5Exception::HDF5Exception(const char *format, ...)
{
    va_list argptr;
    va_start(argptr,format);
    size_t needed = vsnprintf(nullptr, 0, format, argptr) + 1;
    char buf[needed];
    va_end(argptr);

    va_start(argptr,format);
    vsprintf(buf, format, argptr);
    va_end(argptr);

    msg = buf;
}

const char *HDF5Exception::what() const noexcept { return msg.c_str(); }

bool hdf5DatasetExists(hid_t file_id, const std::string &datasetName)
{
    return hdf5DatasetExists(file_id, datasetName);
}

void closeHDF5File(hid_t file_id)
{
    std::lock_guard<mutexHdfType> lock(mutexHdf);
    if(file_id < 1)
        throw HDF5Exception("closeHDF5File: Invalid file handle given.");

    H5_SAVE_ERROR_HANDLER;
    herr_t status = H5Fclose(file_id);

    if (status < 0) {
        error("closeHDF5File failed to close HDF5 file.");
        H5Eprint(H5E_DEFAULT, stderr);
    }
    H5_RESTORE_ERROR_HANDLER;
}

void hdf5EnsureGroupExists(hid_t file_id, std::string const& groupName)
{
    hdf5EnsureGroupExists(file_id, groupName.c_str());
}

void hdf5CreateExtendableString1DArray(hid_t file_id, const char *datasetPath)
{
    int rank = 1;
    hsize_t initialDimensions[1] = {0};
    hsize_t maximumDimensions[1] = {H5S_UNLIMITED};

    std::lock_guard<mutexHdfType> lock(mutexHdf);

    H5::DataSpace dataspace(rank, initialDimensions, maximumDimensions);
    H5::H5File file (file_id);

    // need chunking for extendable dataset
    hsize_t chunkDimensions[1] = {1};
    H5::DSetCreatPropList datasetCreationProperty;
    datasetCreationProperty.setChunk(rank, chunkDimensions);

    H5::StrType strType(0, H5T_VARIABLE);
    RELEASE_ASSERT(H5T_STRING == H5Tget_class(strType.getId())
                   && H5Tis_variable_str(strType.getId()), "");

    auto dataset = file.createDataSet(datasetPath, strType, dataspace,
                                      datasetCreationProperty);
}

void hdf5ExtendAndWriteToString1DArray(hid_t file_id, const char *datasetPath,
                                       const std::string &buffer)
{
    std::lock_guard<mutexHdfType> lock(mutexHdf);

    H5::H5File file (file_id);
    auto dataset = file.openDataSet(datasetPath);

    // extend
    auto filespace = dataset.getSpace();
    int rank = filespace.getSimpleExtentNdims();
    if(rank != 1)
        throw HDF5Exception("Only works for 1D arrays!");

    hsize_t currentDimensions[1];
    filespace.getSimpleExtentDims(currentDimensions);
    hsize_t newDimensions[1] = {currentDimensions[0] + 1};
    dataset.extend(newDimensions);

    filespace = dataset.getSpace();
    hsize_t offset[1] = {currentDimensions[0]};
    hsize_t slabsize[1] = {1};
    filespace.selectHyperslab(H5S_SELECT_SET, slabsize, offset);

    H5::StrType strType(0, H5T_VARIABLE);
    H5::DataSpace memspace(rank, slabsize);

    dataset.write(buffer, strType, memspace, filespace);
}

void hdf5CreateOrExtendAndWriteToString1DArray(hid_t file_id,
                                               const char *parentPath,
                                               const char *datasetName,
                                               const std::string &buffer)
{
    hdf5EnsureGroupExists(file_id, parentPath);

    std::string fullDatasetPath = std::string(parentPath) + "/" + datasetName;

    if (!hdf5DatasetExists(file_id, fullDatasetPath)) {
        hdf5CreateExtendableString1DArray(file_id, fullDatasetPath.c_str());
    }

    hdf5ExtendAndWriteToString1DArray(file_id, fullDatasetPath.c_str(), buffer);
}

H5::H5File hdf5OpenForReading(const std::string &hdf5Filename)
{
    auto lock = hdf5MutexGetLock();

    H5_SAVE_ERROR_HANDLER;
    try {
        auto file = H5::H5File(hdf5Filename, H5F_ACC_RDONLY);
        H5_RESTORE_ERROR_HANDLER;
        return file;
    } catch (...) {
        logmessage(LOGLVL_CRITICAL,
                   "failed to open HDF5 file '%s'.",
                   hdf5Filename.c_str());
        printBacktrace(20);
        H5Ewalk2(H5E_DEFAULT, H5E_WALK_DOWNWARD, hdf5ErrorStackWalker_cb,
                 nullptr);
        H5_RESTORE_ERROR_HANDLER;
        throw(HDF5Exception());
    }
}

void hdf5EnsureGroupExists(const H5::H5File & file, const std::string &groupName)
{
    hdf5EnsureGroupExists(file, groupName);
}

bool hdf5DatasetExists(const H5::H5File &file, const std::string &datasetName)
{
    return hdf5DatasetExists(file, datasetName);
}

} // namespace parpe
