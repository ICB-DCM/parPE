/**
 * Functions for HDF5 I/O
 *
 * NOTE: Use only `const char*` versions of any HDF5 functions, not the
 * `std::string` version. On many systems, HDF5 libraries are still not compiled
 * with C++11 support, but use the old C++03 ABI, which will lead to linking
 * issues.
 */

#include <parpecommon/hdf5Misc.h>
#include <parpecommon/logging.h>
#include <parpecommon/misc.h>

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <utility>
#include <unistd.h>
#include <sys/stat.h>

#ifndef __INTEL_COMPILER
#include <filesystem>
using std::filesystem::path;
using std::filesystem::create_directories;
#else
#include <experimental/filesystem>
using std::experimental::filesystem::path;
using std::experimental::filesystem::create_directories;
#endif

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
    Ensures(err_desc != nullptr);
    const int indent = 2;

    std::unique_ptr<char, decltype(std::free) *>
            maj_str { H5Eget_major(err_desc->maj_num), &std::free };
    std::unique_ptr<char, decltype(std::free) *>
            min_str { H5Eget_minor(err_desc->min_num), &std::free };

    logmessage(LOGLVL_CRITICAL, "%*s#%03d: %s line %u in %s(): %s", indent, "",
               n, err_desc->file_name, err_desc->line, err_desc->func_name,
               err_desc->desc);
    logmessage(LOGLVL_CRITICAL, "%*smajor(%02d): %s", indent * 2, "",
               err_desc->maj_num, maj_str.get());
    logmessage(LOGLVL_CRITICAL, "%*sminor(%02d): %s", indent * 2, "",
               err_desc->min_num, min_str.get());

    return 0;
}


bool hdf5GroupExists(H5::H5File const& file, const std::string &groupName)
{
    std::lock_guard<mutexHdfType> lock(mutexHdf);

    // switch off error handler, check existence and re-enable
    H5_SAVE_ERROR_HANDLER;

    herr_t status = H5Gget_objinfo(file.getId(), groupName.c_str(), false, nullptr);

    H5_RESTORE_ERROR_HANDLER;

    return status >= 0;
}

void hdf5EnsureGroupExists(H5::H5File const& file,
                           const std::string &groupName) {
    if (!hdf5GroupExists(file, groupName)) {
        hdf5CreateGroup(file, groupName, true);
    }
}

void hdf5CreateGroup(const H5::H5File &file, const std::string &groupPath, bool recursively)
{
    std::lock_guard<mutexHdfType> lock(mutexHdf);

    // requires HDF5 >=1.10.6, so needs some C here
    // groupCreationPropertyList.setCreateIntermediateGroup(recursively);
    auto groupCreationPropertyListTmp = H5Pcreate(H5P_LINK_CREATE);
    H5Pset_create_intermediate_group(groupCreationPropertyListTmp, recursively);
    H5::LinkCreatPropList groupCreationPropertyList(groupCreationPropertyListTmp);

    try {
        auto group = file.createGroup(groupPath.c_str(), groupCreationPropertyList);
    }  catch (H5::Exception const&) {
        throw(HDF5Exception("Failed to create group in hdf5CreateGroup:" +
                            groupPath));
    }
}

void hdf5CreateExtendableDouble2DArray(const H5::H5File &file,
                                       const std::string &datasetPath,
                                       hsize_t stride)
{
    constexpr int rank = 2;
    hsize_t initialDimensions[2] = {stride, 0};
    hsize_t maximumDimensions[2] = {stride, H5S_UNLIMITED};

    std::lock_guard<mutexHdfType> lock(mutexHdf);

    H5::DataSpace dataspace(rank, initialDimensions, maximumDimensions);

    // need chunking for extendable dataset
    hsize_t chunkDimensions[2] = {stride, 1};
    H5::DSetCreatPropList dSetCreatPropList;
    dSetCreatPropList.setChunk(rank, chunkDimensions);

    auto dataset = file.createDataSet(datasetPath.c_str(), H5::PredType::NATIVE_DOUBLE,
                                      dataspace, dSetCreatPropList);
}

void hdf5Extend2ndDimensionAndWriteToDouble2DArray(const H5::H5File &file, const std::string &datasetPath, gsl::span<const double> buffer)
{
    std::lock_guard<mutexHdfType> lock(mutexHdf);

    auto dataset = file.openDataSet(datasetPath.c_str());

    // check rank
    auto filespace = dataset.getSpace();
    int rank = filespace.getSimpleExtentNdims();
    if (rank != 2) {
        throw HDF5Exception("Failed to write data in "
                            "hdf5Extend2ndDimensionAndWriteToDouble2DArray: "
                            "not of rank 2 (%d) when writing %s",
                            rank, datasetPath.c_str());
    }

    // extend
    hsize_t currentDimensions[2];
    filespace.getSimpleExtentDims(currentDimensions);

    Expects(buffer.size() == currentDimensions[0]);

    hsize_t newDimensions[2] = {currentDimensions[0],
                                currentDimensions[1] + 1};
    dataset.extend(newDimensions);

    filespace = dataset.getSpace();
    hsize_t offset[2] = {0, currentDimensions[1]};
    hsize_t slabsize[2] = {currentDimensions[0], 1};

    filespace.selectHyperslab(H5S_SELECT_SET, slabsize, offset);

    H5::DataSpace memspace(rank, slabsize);
    dataset.write(buffer.data(), H5::PredType::NATIVE_DOUBLE, memspace, filespace);
}

void hdf5Extend3rdDimensionAndWriteToDouble3DArray(const H5::H5File &file,
                                                   std::string const& datasetPath,
                                                   gsl::span<const double> buffer)
{
    std::lock_guard<mutexHdfType> lock(mutexHdf);

    auto dataset = file.openDataSet(datasetPath.c_str());

    // extend
    auto filespace = dataset.getSpace();
    int rank = filespace.getSimpleExtentNdims();
    RELEASE_ASSERT(rank == 3, "Only works for 3D arrays!");

    hsize_t currentDimensions[3];
    filespace.getSimpleExtentDims(currentDimensions);

    hsize_t newDimensions[3] = {currentDimensions[0],
                                currentDimensions[1],
                                currentDimensions[2] + 1};
    dataset.extend(newDimensions);

    filespace = dataset.getSpace();
    hsize_t offset[3] = {0, 0, currentDimensions[2]};
    hsize_t slabsize[3] = {currentDimensions[0], currentDimensions[1], 1};

    filespace.selectHyperslab(H5S_SELECT_SET, slabsize, offset);

    auto memspace = H5::DataSpace(rank, slabsize);

    dataset.write(buffer.data(), H5::PredType::NATIVE_DOUBLE, memspace, filespace);
}

void hdf5CreateOrExtendAndWriteToDouble2DArray(const H5::H5File &file,
                                               const std::string &parentPath,
                                               const std::string &datasetName,
                                               gsl::span<const double> buffer)
{
    hdf5EnsureGroupExists(file, parentPath);

    std::string fullDatasetPath = std::string(parentPath) + "/" + datasetName;

    if (!file.nameExists(fullDatasetPath.c_str())) {
        hdf5CreateExtendableDouble2DArray(
                    file, fullDatasetPath, buffer.size());
    }

    hdf5Extend2ndDimensionAndWriteToDouble2DArray(
                file, fullDatasetPath, buffer);
}

void hdf5CreateOrExtendAndWriteToDouble3DArray(const H5::H5File &file,
                                               const std::string &parentPath,
                                               const std::string &datasetName,
                                               gsl::span<const double> buffer,
                                               hsize_t stride1,
                                               hsize_t stride2) {
    hdf5EnsureGroupExists(file, parentPath);

    std::string fullDatasetPath = std::string(parentPath) + "/" + datasetName;

    if (!file.nameExists(fullDatasetPath.c_str())) {
        hdf5CreateExtendableDouble3DArray(
                    file, fullDatasetPath, stride1, stride2);
    }

    hdf5Extend3rdDimensionAndWriteToDouble3DArray(file,
                                                  fullDatasetPath,
                                                  buffer);

}

void hdf5CreateOrExtendAndWriteToInt2DArray(const H5::H5File &file,
                                            const std::string &parentPath,
                                            const std::string &datasetName,
                                            gsl::span<const int> buffer)
{
    std::lock_guard<mutexHdfType> lock(mutexHdf);

    hdf5EnsureGroupExists(file, parentPath);

    auto fullDatasetPath = std::string(parentPath) + "/" + datasetName;

    if (!file.nameExists(fullDatasetPath.c_str())) {
        hdf5CreateExtendableInt2DArray(
                    file, fullDatasetPath, buffer.size());
    }

    hdf5Extend2ndDimensionAndWriteToInt2DArray(file, fullDatasetPath,
                                               buffer);
}

void hdf5Extend2ndDimensionAndWriteToInt2DArray(const H5::H5File &file,
                                                const std::string &datasetPath,
                                                gsl::span<const int> buffer)
{
    std::lock_guard<mutexHdfType> lock(mutexHdf);

    auto dataset = file.openDataSet(datasetPath.c_str());

    // extend
    auto filespace = dataset.getSpace();
    int rank = filespace.getSimpleExtentNdims();
    if(rank != 2)
        throw HDF5Exception("Only works for 2D arrays!");

    hsize_t currentDimensions[2];
    filespace.getSimpleExtentDims(currentDimensions);
    Expects(buffer.size() == currentDimensions[0]);

    hsize_t newDimensions[2] = {currentDimensions[0],
                                currentDimensions[1] + 1};
    dataset.extend(newDimensions);

    filespace = dataset.getSpace();
    hsize_t offset[2] = {0, currentDimensions[1]};
    hsize_t slabsize[2] = {currentDimensions[0], 1};

    filespace.selectHyperslab(H5S_SELECT_SET, slabsize, offset);

    H5::DataSpace memspace(rank, slabsize);
    dataset.write(buffer.data(), H5::PredType::NATIVE_INT, memspace, filespace);
}

void hdf5CreateExtendableInt2DArray(const H5::H5File &file,
                                    const std::string &datasetPath,
                                    hsize_t stride)
{
    std::lock_guard<mutexHdfType> lock(mutexHdf);

    int rank = 2;
    hsize_t initialDimensions[2] = {stride, 0};
    hsize_t maximumDimensions[2] = {stride, H5S_UNLIMITED};

    H5::DataSpace dataspace(rank, initialDimensions, maximumDimensions);

    // need chunking for extendable dataset
    hsize_t chunkDimensions[2] = {stride, 1};
    H5::DSetCreatPropList datasetCreationProperty;
    datasetCreationProperty.setChunk(rank, chunkDimensions);

    Expects(H5Tget_size(H5T_NATIVE_INT) == sizeof(int));
    auto dataset = file.createDataSet(datasetPath.c_str(), H5::PredType::NATIVE_INT,
                                      dataspace, datasetCreationProperty);

}

void hdf5CreateExtendableDouble3DArray(const H5::H5File &file,
                                       const std::string &datasetPath,
                                       hsize_t stride1,
                                       hsize_t stride2)
{

    int rank = 3;
    hsize_t initialDimensions[3] = {stride1, stride2, 0};
    hsize_t maximumDimensions[3] = {stride1, stride2, H5S_UNLIMITED};

    std::lock_guard<mutexHdfType> lock(mutexHdf);

    H5::DataSpace dataspace(rank, initialDimensions, maximumDimensions);

    // need chunking for extendable dataset
    hsize_t chunkDimensions[3] = {stride1, stride2, 1};
    H5::DSetCreatPropList datasetCreationProperty;
    datasetCreationProperty.setChunk(rank, chunkDimensions);

    auto dataset = file.createDataSet(datasetPath.c_str(), H5::PredType::NATIVE_DOUBLE,
                                      dataspace, datasetCreationProperty);
}

void hdf5Read2DDoubleHyperslab(const H5::H5File &file,
                               const std::string &path,
                               hsize_t size0,
                               hsize_t size1,
                               hsize_t offset0,
                               hsize_t offset1,
                               gsl::span<double> buffer)
{
    Expects(buffer.size() == size0 * size1);

    std::lock_guard<mutexHdfType> lock(mutexHdf);

    auto dataset = file.openDataSet(path.c_str());
    auto dataspace = dataset.getSpace();
    hsize_t offset[] = {offset0, offset1};
    hsize_t count[] = {size0, size1};

    const int ndims = dataspace.getSimpleExtentNdims();
    RELEASE_ASSERT(ndims == 2 && "Only works for 2D arrays!", "");
    hsize_t dims[ndims];
    dataspace.getSimpleExtentDims(dims);
    // printf("%lld %lld, %lld %lld, %lld %lld\n", dims[0], dims[1], offset0,
    // offset1, size0, size1);
    RELEASE_ASSERT(dims[0] >= offset0 && dims[0] >= size0,
            "Offset larger than dataspace dimensions!");
    RELEASE_ASSERT(dims[1] >= offset1 && dims[1] >= size1,
            "Offset larger than dataspace dimensions!");

    dataspace.selectHyperslab(H5S_SELECT_SET, count, offset);

    H5::DataSpace memspace(2, count);

    dataset.read(buffer.data(), H5::PredType::NATIVE_DOUBLE, memspace, dataspace);
}

std::vector<int> hdf5Read1DIntegerHyperslab(H5::H5File const& file,
                                            std::string const& path,
                                            hsize_t count,
                                            hsize_t offset)
{
    std::lock_guard<mutexHdfType> lock(mutexHdf);

    H5::DataSet dataset = file.openDataSet(path.c_str());
    H5::DataSpace filespace = dataset.getSpace();

    const int ndims = filespace.getSimpleExtentNdims();
    RELEASE_ASSERT(ndims == 1, "Only works for 1D arrays!");
    hsize_t length;
    filespace.getSimpleExtentDims(&length);

    RELEASE_ASSERT(length >= offset,
                   "Offset larger than dataspace dimensions!");

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

    H5::DataSet dataset = file.openDataSet(path.c_str());
    H5::DataSpace filespace = dataset.getSpace();

    const int ndims = filespace.getSimpleExtentNdims();
    RELEASE_ASSERT(ndims == 2, "Only works for 2D arrays!");
    hsize_t dims[ndims];
    filespace.getSimpleExtentDims(dims);

    hsize_t offset[] = {offset0, offset1};
    hsize_t count[] = {size0, size1};
    // printf("%lld %lld, %lld %lld, %lld %lld\n", dims[0], dims[1], offset0,
    // offset1, size0, size1);
    if(offset0 >= dims[0] || size0 > dims[0] || offset1 >= dims[1]
            || size1 > dims[1]) {
        std::stringstream ss;
        ss << "Offset larger than dataspace dimensions! " << "dims: " << dims[0]
           << "," << dims[1] << " offsets: " << offset0 << "," << offset1
           << " size: " << size0 << "," << size1;
        printBacktrace();
        throw HDF5Exception(ss.str());
    }

    filespace.selectHyperslab(H5S_SELECT_SET, count, offset);

    H5::DataSpace memspace(2, count);
    std::vector<int> buffer(size0 * size1);

    dataset.read(buffer.data(), H5::PredType::NATIVE_INT, memspace, filespace);

    return buffer;
}

void hdf5Read3DDoubleHyperslab(H5::H5File const& file,
                               std::string const& path,
                               hsize_t size0,
                               hsize_t size1,
                               hsize_t size2,
                               hsize_t offset0,
                               hsize_t offset1,
                               hsize_t offset2,
                               gsl::span<double> buffer)
{
    std::lock_guard<mutexHdfType> lock(mutexHdf);

    const int rank = 3;
    auto dataset = file.openDataSet(path.c_str());
    auto dataspace = dataset.getSpace();
    hsize_t offset[] = {offset0, offset1, offset2};
    hsize_t count[] = {size0, size1, size2};

    const int ndims = dataspace.getSimpleExtentNdims();
    RELEASE_ASSERT(ndims == rank, "Only works for 3D arrays!");
    hsize_t dims[ndims];
    dataspace.getSimpleExtentDims(dims);
    RELEASE_ASSERT(dims[0] >= offset0 && dims[0] >= size0,
            "Offset larger than dataspace dimensions!");
    RELEASE_ASSERT(dims[1] >= offset1 && dims[1] >= size1,
            "Offset larger than dataspace dimensions!");
    RELEASE_ASSERT(dims[2] >= offset2 && dims[2] >= size2,
            "Offset larger than dataspace dimensions!");

    dataspace.selectHyperslab(H5S_SELECT_SET, count, offset);

    H5::DataSpace memspace(rank, count);
    dataset.read(buffer.data(), H5::PredType::NATIVE_DOUBLE, memspace, dataspace);
}


std::vector<double> hdf5Get3DDoubleHyperslab(const H5::H5File &file, const std::string &path,
                                             hsize_t size0, hsize_t size1,
                                             hsize_t size2,
                                             hsize_t offset0, hsize_t offset1,
                                             hsize_t offset2)
{
    std::vector<double> buffer(size0 * size1 * size2);
    hdf5Read3DDoubleHyperslab(file, path, size0, size1, size2,
                              offset0, offset1, offset2, buffer);
    return buffer;
}


bool hdf5AttributeExists(const H5::H5File &file,
                         const std::string &datasetPath,
                         const std::string &attributeName)
{
    std::lock_guard<mutexHdfType> lock(mutexHdf);

    int exists = false;

    H5_SAVE_ERROR_HANDLER;

    auto loc = H5Oopen(file.getId(), datasetPath.c_str(), H5P_DEFAULT);
    if (loc >= 0) {
        exists = H5LTfind_attribute(loc, attributeName.c_str());
        H5Oclose(loc);
    }

    H5_RESTORE_ERROR_HANDLER;

    return exists;
}

void hdf5WriteStringAttribute(const H5::H5File &file,
                              const std::string &datasetPath,
                              const std::string &attributeName,
                              const std::string &attributeValue)
{
    std::lock_guard<mutexHdfType> lock(mutexHdf);

    int ret = H5LTset_attribute_string(
                file.getId(), datasetPath.c_str(),
                attributeName.c_str(), attributeValue.c_str());
    if(ret < 0)
        throw HDF5Exception("Unable to write attribute %s on %s",
                            datasetPath.c_str(), attributeName.c_str());
}

H5::H5File hdf5CreateFile(const std::string &filename,
                          bool overwrite)
{
    // Create parent folders
    path dirname(filename);
    dirname.remove_filename();
    if(!dirname.empty())
        create_directories(dirname);

    std::lock_guard<mutexHdfType> lock(mutexHdf);

    if (!overwrite) {
        struct stat st = {};
        bool fileExists = stat(filename.c_str(), &st) == 0;

        if(fileExists)
            throw HDF5Exception("Result file exists " + filename);
    }

    try {
        return H5::H5File(filename.c_str(), H5F_ACC_TRUNC);
    }  catch (H5::Exception const& e) {
        printBacktrace();
        throw HDF5Exception("hdf5CreateFile: Failed to create file %s. "
                            "Is this file opened by another process?",
                            filename.c_str());
    }
}


void hdf5GetDatasetDimensions(const H5::H5File &file,
                              const std::string &path,
                              hsize_t nDimsExpected,
                              int *d1,
                              int *d2,
                              int *d3,
                              int *d4)
{
    std::lock_guard<mutexHdfType> lock(mutexHdf);
    H5_SAVE_ERROR_HANDLER;

    auto dataset = file.openDataSet(path.c_str());
    auto dataspace = dataset.getSpace();

    const int nDimsActual = dataspace.getSimpleExtentNdims();
    if(nDimsActual != (signed)nDimsExpected)
        throw(HDF5Exception("Dataset rank (" + std::to_string(nDimsActual)
                            + ") does not match nDims argument ("
                            + std::to_string(nDimsExpected) + ")"));

    hsize_t dims[nDimsExpected];
    dataspace.getSimpleExtentDims(dims);

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



void hdf5CreateExtendableString1DArray(const H5::H5File &file, const std::string &datasetPath)
{
    int rank = 1;
    hsize_t initialDimensions[1] = {0};
    hsize_t maximumDimensions[1] = {H5S_UNLIMITED};

    std::lock_guard<mutexHdfType> lock(mutexHdf);

    H5::DataSpace dataspace(rank, initialDimensions, maximumDimensions);

    // need chunking for extendable dataset
    hsize_t chunkDimensions[1] = {1};
    H5::DSetCreatPropList datasetCreationProperty;
    datasetCreationProperty.setChunk(rank, chunkDimensions);

    H5::StrType strType(0, H5T_VARIABLE);
    Expects(H5T_STRING == H5Tget_class(strType.getId())
            && H5Tis_variable_str(strType.getId()));

    file.createDataSet(datasetPath.c_str(), strType,
                       dataspace, datasetCreationProperty);
}

void hdf5ExtendAndWriteToString1DArray(const H5::H5File &file, const std::string &datasetPath,
                                       const std::string &buffer)
{
    std::lock_guard<mutexHdfType> lock(mutexHdf);

    auto dataset = file.openDataSet(datasetPath.c_str());

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

    dataset.write(buffer.c_str(), strType, memspace, filespace);
}

void hdf5CreateOrExtendAndWriteToString1DArray(const H5::H5File &file,
                                               const std::string &parentPath,
                                               const std::string &datasetName,
                                               const std::string &buffer)
{
    hdf5EnsureGroupExists(file, parentPath);

    std::string fullDatasetPath = std::string(parentPath) + "/" + datasetName;

    if (!file.nameExists(fullDatasetPath.c_str())) {
        hdf5CreateExtendableString1DArray(file, fullDatasetPath);
    }

    hdf5ExtendAndWriteToString1DArray(file, fullDatasetPath, buffer);
}

H5::H5File hdf5OpenForReading(const std::string &hdf5Filename)
{
    [[maybe_unused]] auto lock = hdf5MutexGetLock();

    H5_SAVE_ERROR_HANDLER;
    try {
        auto file = H5::H5File(hdf5Filename.c_str(), H5F_ACC_RDONLY);
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
        throw(HDF5Exception("Unable to open HDF5 file " + hdf5Filename));
    }
}

H5::H5File hdf5OpenForAppending(const std::string &hdf5Filename)
{
    [[maybe_unused]] auto lock = hdf5MutexGetLock();

    H5::H5File file;
    H5_SAVE_ERROR_HANDLER;

    // Open for append or create
    try {
        file = H5::H5File(hdf5Filename.c_str(), H5F_ACC_RDWR);
        H5_RESTORE_ERROR_HANDLER;
    } catch (H5::FileIException const&) {
        H5_RESTORE_ERROR_HANDLER;
        // create if doesn't exist
        file = H5::H5File(hdf5Filename.c_str(), H5F_ACC_EXCL);
    }

    return file;
}

std::vector<std::string> hdf5Read1dStringDataset(
        H5::H5File const& file, const std::string &datasetPath)
{
    [[maybe_unused]] auto lock = hdf5MutexGetLock();
    auto dataset = file.openDataSet(datasetPath.c_str());
    auto filespace = dataset.getSpace();

    const int ndims = filespace.getSimpleExtentNdims();
    RELEASE_ASSERT(ndims == 1, "Only works for 1D arrays!");

    auto dtype = dataset.getDataType();
    auto native_type = H5Tget_native_type(dtype.getId(), H5T_DIR_DEFAULT);
    H5::StrType tid1(0, H5T_VARIABLE);
    if(!H5Tequal(native_type, tid1.getId()))
        throw HDF5Exception("Data type mismatch");

    hsize_t length;
    filespace.getSimpleExtentDims(&length);
    std::vector<char*> buffer(length);
    dataset.read((void*)buffer.data(), dtype);

    std::vector<std::string> strBuffer(buffer.size());
    for(int i = 0; i < (int) buffer.size(); ++i) {
        strBuffer[i] = buffer[i];
    }
    return strBuffer;
}

void hdf5Write1dStringDataset(
        const H5::H5File &file, const std::string &parentPath,
        const std::string &datasetPath, std::vector<std::string> const& buffer)
{
    [[maybe_unused]] auto lock = hdf5MutexGetLock();

    const int dims = 1;
    hsize_t dims0 = buffer.size();
    H5::DataSpace sid1(dims, &dims0);

    H5::StrType tid1(0, H5T_VARIABLE);
    RELEASE_ASSERT(H5T_STRING == H5Tget_class(tid1.getId())
                   || !H5Tis_variable_str(tid1.getId()), "String type failure.");

    hdf5EnsureGroupExists(file, parentPath);
    std::string fullpath(parentPath + "/" + datasetPath);
    auto dataset = file.createDataSet(fullpath.c_str(), tid1, sid1);

    // we need character pointers
    std::vector<const char *> charPtrBuffer(buffer.size());
    for(int i = 0; i < (int) buffer.size(); ++i) {
        charPtrBuffer[i] = buffer[i].c_str();
    }
    dataset.write((void*)charPtrBuffer.data(), tid1);
}

} // namespace parpe
