#ifndef HDF5_MISC_H
#define HDF5_MISC_H

#include <hdf5.h>
#include <hdf5_hl.h>

#include <pthread.h>
#include <stdbool.h>
#include <exception>
#include <string>

namespace parPE {

class HDF5Exception : public std::exception {
public:
    HDF5Exception(std::string msg = "") : msg(msg) {}

    const char* what() const noexcept override { return msg.c_str(); }

    std::string msg;

};


void initHDF5Mutex();

void hdf5LockMutex();

void hdf5UnlockMutex();

void destroyHDF5Mutex();

#define H5_SAVE_ERROR_HANDLER                                                  \
    herr_t (*old_func)(void *);                                                \
    void *old_client_data;                                                     \
    H5Eget_auto1(&old_func, &old_client_data);                                 \
    H5Eset_auto1(NULL, NULL)

#define H5_RESTORE_ERROR_HANDLER H5Eset_auto1(old_func, old_client_data)

herr_t
hdf5ErrorStackWalker_cb(unsigned int n, const H5E_error_t *err_desc,
                        void *client_data); // TODO: also use for resultwriter

bool hdf5DatasetExists(hid_t file_id, const char *datasetName);

bool hdf5GroupExists(hid_t file_id, const char *groupName);

bool hdf5EnsureGroupExists(hid_t file_id, const char *groupName);

void hdf5CreateGroup(hid_t file_id, const char *groupPath, bool recursively);

hid_t hdf5OpenFile(const char *filename,
                   bool overwrite);

void hdf5CreateExtendableDouble2DArray(hid_t file_id, const char *datasetPath,
                                       hsize_t stride);

void hdf5CreateExtendableInt2DArray(hid_t file_id, const char *datasetPath,
                                    hsize_t stride);

void hdf5CreateExtendableDouble3DArray(hid_t file_id, const char *datasetPath,
                                       hsize_t stride1, hsize_t stride2);

void hdf5Extend2ndDimensionAndWriteToDouble2DArray(hid_t file_id,
                                                   const char *datasetPath,
                                                   const double *buffer);

void hdf5Extend2ndDimensionAndWriteToInt2DArray(hid_t file_id,
                                                const char *datasetPath,
                                                const int *buffer);

void hdf5CreateOrExtendAndWriteToDouble2DArray(hid_t file_id,
                                               const char *parentPath,
                                               const char *datasetName,
                                               const double *buffer,
                                               hsize_t stride);

void hdf5CreateOrExtendAndWriteToInt2DArray(hid_t file_id,
                                            const char *parentPath,
                                            const char *datasetName,
                                            const int *buffer, hsize_t stride);

void hdf5CreateOrExtendAndWriteToDouble3DArray(hid_t file_id,
                                               const char *parentPath,
                                               const char *datasetName,
                                               const double *buffer,
                                               hsize_t stride1, hsize_t stride2);

int hdf5Read2DDoubleHyperslab(hid_t file_id, const char *path, hsize_t size0,
                              hsize_t size1, hsize_t offset0, hsize_t offset1,
                              double *buffer);

int hdf5Read3DDoubleHyperslab(hid_t file_id, const char *path, hsize_t size0,
                              hsize_t size1, hsize_t size2, hsize_t offset0,
                              hsize_t offset1, hsize_t offset2, double *buffer);

void hdf5GetDatasetDimensions2D(hid_t file_id, const char *path, int *d1,
                                int *d2);

void hdf5GetDatasetDimensions3D(hid_t file_id, const char *path, int *d1,
                                int *d2, int *d3);

int hdf5AttributeExists(hid_t fileId, const char *datasetPath,
                        const char *attributeName);

int hdf5WriteStringAttribute(hid_t fileId, const char *datasetPath,
                             const char *attributeName,
                             const char *attributeValue);

} // namespace parPE
#endif
