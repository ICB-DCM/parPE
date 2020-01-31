#ifndef HDF5_MISC_H
#define HDF5_MISC_H
#include <parpecommon/misc.h>

#include <hdf5.h>
#include <hdf5_hl.h>
#include <H5Cpp.h>

#include <pthread.h>
#include <exception>
#include <string>
#include <mutex>
#include <cstdarg>

namespace parpe {

class HDF5Exception : public std::exception {
public:
    HDF5Exception(std::string msg = "");

    HDF5Exception(const char *format, ...);

    const char* what() const noexcept;

    std::string msg;
    std::string stackTrace;
};


typedef std::recursive_mutex mutexHdfType;

void initHDF5Mutex();

std::unique_lock<mutexHdfType> hdf5MutexGetLock();

#define H5_SAVE_ERROR_HANDLER                                                  \
    herr_t (*old_func)(void *);                                                \
    void *old_client_data;                                                     \
    H5Eget_auto1(&old_func, &old_client_data);                                 \
    H5Eset_auto1(nullptr, nullptr)

#define H5_RESTORE_ERROR_HANDLER H5Eset_auto1(old_func, old_client_data)

herr_t
hdf5ErrorStackWalker_cb(unsigned int n, const H5E_error_t *err_desc, void *); // TODO: also use for resultwriter

bool hdf5DatasetExists(hid_t file_id, std::string const& datasetName);

bool hdf5DatasetExists(hid_t file_id, const char *datasetName);

bool hdf5DatasetExists(H5::H5File const& file, const std::string &datasetName);

bool hdf5GroupExists(hid_t file_id, const char *groupName);

void hdf5EnsureGroupExists(hid_t file_id, const char *groupName);

void hdf5EnsureGroupExists(hid_t file_id, const std::string &groupName);

void hdf5EnsureGroupExists(H5::H5File const& file,
                           const std::string &groupName);

void hdf5CreateGroup(hid_t file_id, const char *groupPath, bool recursively = false);

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
hid_t hdf5CreateFile(const char *filename,
                   bool overwrite = false);

H5::H5File hdf5OpenForReading(std::string const& hdf5Filename);

void closeHDF5File(hid_t file_id);

void hdf5CreateExtendableDouble2DArray(hid_t file_id, const char *datasetPath,
                                       hsize_t stride);

void hdf5CreateExtendableInt2DArray(hid_t file_id, const char *datasetPath,
                                    hsize_t stride);

void hdf5CreateExtendableDouble3DArray(hid_t file_id, const char *datasetPath,
                                       hsize_t stride1, hsize_t stride2);

void hdf5CreateExtendableString1DArray(hid_t file_id, const char *datasetPath);

void hdf5Extend2ndDimensionAndWriteToDouble2DArray(hid_t file_id,
                                                   const char *datasetPath,
                                                   gsl::span<const double> buffer);

void hdf5Extend2ndDimensionAndWriteToInt2DArray(hid_t file_id,
                                                const char *datasetPath,
                                                gsl::span<const int> buffer);

void hdf5ExtendAndWriteToString1DArray(hid_t file_id,
                                       const char *datasetPath,
                                       std::string const& buffer);

void hdf5CreateOrExtendAndWriteToDouble2DArray(hid_t file_id,
                                               const char *parentPath,
                                               const char *datasetName,
                                               gsl::span<const double> buffer);

void hdf5CreateOrExtendAndWriteToInt2DArray(hid_t file_id,
                                            const char *parentPath,
                                            const char *datasetName,
                                            gsl::span<const int> buffer);

void hdf5CreateOrExtendAndWriteToDouble3DArray(hid_t file_id,
                                               const char *parentPath,
                                               const char *datasetName,
                                               gsl::span<const double> buffer,
                                               hsize_t stride1,
                                               hsize_t stride2);

void hdf5CreateOrExtendAndWriteToString1DArray(hid_t file_id,
                                               const char *parentPath,
                                               const char *datasetName,
                                               std::string const& buffer);

int hdf5Read2DDoubleHyperslab(hid_t file_id, const char *path, hsize_t size0,
                              hsize_t size1, hsize_t offset0, hsize_t offset1,
                              gsl::span<double> buffer);

int hdf5Read3DDoubleHyperslab(hid_t file_id, const char *path, hsize_t size0,
                              hsize_t size1, hsize_t size2, hsize_t offset0,
                              hsize_t offset1, hsize_t offset2,
                              gsl::span<double> buffer);

std::vector<double> hdf5Get3DDoubleHyperslab(hid_t file_id, const char *path, hsize_t size0,
                              hsize_t size1, hsize_t size2, hsize_t offset0,
                              hsize_t offset1, hsize_t offset2);

std::vector<int> hdf5Read1DIntegerHyperslab(
        const H5::H5File &file, std::string const& path,
        hsize_t count, hsize_t offset);

std::vector<int> hdf5Read2DIntegerHyperslab(
        H5::H5File const& file, std::string const& path,
        hsize_t size0, hsize_t size1, hsize_t offset0, hsize_t offset1);

void hdf5GetDatasetDimensions(hid_t file_id, const char *path,
                              hsize_t nDimsExpected,
                              int *d1 = nullptr, int *d2 = nullptr,
                              int *d3 = nullptr, int *d4 = nullptr);

bool hdf5AttributeExists(hid_t fileId, const char *datasetPath,
                        const char *attributeName);

void hdf5WriteStringAttribute(hid_t fileId, const char *datasetPath,
                             const char *attributeName,
                             const char *attributeValue);

std::vector<std::string> hdf5Read1dStringDataset(
        const H5::H5File &file, std::string const& datasetPath);

} // namespace parpe
#endif
