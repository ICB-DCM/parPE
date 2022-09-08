#ifndef HDF5_MISC_H
#define HDF5_MISC_H
#include <parpecommon/misc.h>

#include <hdf5.h>
#include <hdf5_hl.h>
#include <H5Cpp.h>

#include <exception>
#include <string>
#include <mutex>
#include <cstdarg>

namespace parpe {

class HDF5Exception : public std::exception {
public:
    explicit HDF5Exception(std::string msg = "");

    explicit HDF5Exception(const char *format, ...);

    const char* what() const noexcept override;

    std::string msg;
    std::string stackTrace;
};


using mutexHdfType = std::recursive_mutex;

void initHDF5Mutex();

std::unique_lock<mutexHdfType> hdf5MutexGetLock();

#define H5_SAVE_ERROR_HANDLER                                                  \
    herr_t (*old_func)(void *);                                                \
    void *old_client_data;                                                     \
    H5Eget_auto1(&old_func, &old_client_data);                                 \
    H5Eset_auto1(nullptr, nullptr)

#define H5_RESTORE_ERROR_HANDLER H5Eset_auto1(old_func, old_client_data)

herr_t
hdf5ErrorStackWalker_cb(unsigned int n, const H5E_error_t *err_desc, void *);

bool hdf5GroupExists(H5::H5File const& file,
                     const std::string &groupName);

void hdf5EnsureGroupExists(H5::H5File const& file,
                           const std::string &groupName);

void hdf5CreateGroup(H5::H5File const& file, std::string const& groupPath, bool recursively = false);

/**
 * @brief Create and open HDF5 file for writing.
 *
 * Creates parent path if it doesn't exist before creating the file.
 * Throws HDF5Exception on failure.
 * @param filename Filename, optionally including path of the new file
 * @param overwrite Overwrite file if exists. If false and file exists,
 * throws HDF5Exception on failure.
 * @return HDF5 file handle of the created/opened file
 */
H5::H5File hdf5CreateFile(std::string const& filename,
                          bool overwrite = false);

H5::H5File hdf5OpenForReading(std::string const& hdf5Filename);

H5::H5File hdf5OpenForAppending(const std::string &hdf5Filename);

void hdf5CreateExtendableDouble2DArray(H5::H5File const& file, std::string const& datasetPath,
                                       hsize_t stride);

void hdf5CreateExtendableInt2DArray(H5::H5File const& file, std::string const& datasetPath,
                                    hsize_t stride);

void hdf5CreateExtendableDouble3DArray(H5::H5File const& file, std::string const& datasetPath,
                                       hsize_t stride1, hsize_t stride2);

void hdf5CreateExtendableString1DArray(H5::H5File const& file, std::string const& datasetPath);

void hdf5Extend2ndDimensionAndWriteToDouble2DArray(H5::H5File const& file,
                                                   std::string const& datasetPath,
                                                   gsl::span<const double> buffer);

void hdf5Extend2ndDimensionAndWriteToInt2DArray(H5::H5File const& file,
                                                std::string const& datasetPath,
                                                gsl::span<const int> buffer);

void hdf5ExtendAndWriteToString1DArray(H5::H5File const& file,
                                       std::string const& datasetPath,
                                       std::string const& buffer);

void hdf5CreateOrExtendAndWriteToDouble2DArray(H5::H5File const& file,
                                               std::string const& parentPath,
                                               std::string const& datasetName,
                                               gsl::span<const double> buffer);

void hdf5CreateOrExtendAndWriteToInt2DArray(H5::H5File const& file,
                                            std::string const& parentPath,
                                            std::string const& datasetName,
                                            gsl::span<const int> buffer);

void hdf5CreateOrExtendAndWriteToDouble3DArray(H5::H5File const& file,
                                               std::string const& parentPath,
                                               std::string const& datasetName,
                                               gsl::span<const double> buffer,
                                               hsize_t stride1,
                                               hsize_t stride2);
void hdf5Extend3rdDimensionAndWriteToDouble3DArray(const H5::H5File &file,
                                                   std::string const& datasetPath,
                                                   gsl::span<const double> buffer);

void hdf5CreateOrExtendAndWriteToString1DArray(H5::H5File const& file,
                                               std::string const& parentPath,
                                               std::string const& datasetName,
                                               std::string const& buffer);

void hdf5Read2DDoubleHyperslab(H5::H5File const& file, std::string const& path, hsize_t size0,
                              hsize_t size1, hsize_t offset0, hsize_t offset1,
                              gsl::span<double> buffer);

void hdf5Read3DDoubleHyperslab(H5::H5File const& file, std::string const& path, hsize_t size0,
                              hsize_t size1, hsize_t size2, hsize_t offset0,
                              hsize_t offset1, hsize_t offset2,
                              gsl::span<double> buffer);

std::vector<double> hdf5Get3DDoubleHyperslab(H5::H5File const& file, std::string const& path, hsize_t size0,
                                             hsize_t size1, hsize_t size2, hsize_t offset0,
                                             hsize_t offset1, hsize_t offset2);

std::vector<int> hdf5Read1DIntegerHyperslab(
        const H5::H5File &file, std::string const& path,
        hsize_t count, hsize_t offset);

std::vector<int> hdf5Read2DIntegerHyperslab(
        H5::H5File const& file, std::string const& path,
        hsize_t size0, hsize_t size1, hsize_t offset0, hsize_t offset1);

void hdf5GetDatasetDimensions(H5::H5File const& file, std::string const& path,
                              hsize_t nDimsExpected,
                              int *d1 = nullptr, int *d2 = nullptr,
                              int *d3 = nullptr, int *d4 = nullptr);

bool hdf5AttributeExists(H5::H5File const& file, std::string const& datasetPath,
                         std::string const& attributeName);

void hdf5WriteStringAttribute(H5::H5File const& file, std::string const& datasetPath,
                              std::string const& attributeName,
                              std::string const& attributeValue);

std::vector<std::string> hdf5Read1dStringDataset(
        const H5::H5File &file, std::string const& datasetPath);

void hdf5Write1dStringDataset(
        const H5::H5File &file, std::string const& parentPath,
        std::string const& datasetPath, std::vector<std::string> const& buffer);

} // namespace parpe
#endif
