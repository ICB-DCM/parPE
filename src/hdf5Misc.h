#ifndef HDF5_MISC_H
#define HDF5_MISC_H

#include <hdf5.h>
#include <hdf5_hl.h>

#include <pthread.h>
#include <stdbool.h>

void initHDF5Mutex();

void hdf5LockMutex();

void hdf5UnlockMutex();

void destroyHDF5Mutex();

#define H5_SAVE_ERROR_HANDLER                  \
    herr_t (*old_func)(void*);                 \
    void *old_client_data;                     \
    H5Eget_auto1(&old_func, &old_client_data); \
    H5Eset_auto1(NULL, NULL)

#define H5_RESTORE_ERROR_HANDLER H5Eset_auto1(old_func, old_client_data)

herr_t hdf5ErrorStackWalker_cb(unsigned int n, const H5E_error_t *err_desc, void *client_data); // TODO: also use for resultwriter


// malloc version of ami_hdf5.cpp
void getDoubleArrayAttributeC(hid_t file_id, const char* optionsObject, const char* attributeName, double **destination, hsize_t *length);
void getIntArrayAttributeC(hid_t file_id, const char* optionsObject, const char* attributeName, int **destination, hsize_t *length);


bool hdf5DatasetExists(hid_t file_id, const char *datasetName);

bool hdf5GroupExists(hid_t file_id, const char *groupName);

void hdf5CreateGroup(hid_t file_id, const char *groupPath, bool recursively);

void hdf5CreateExtendableDouble2DArray(hid_t file_id, const char* datasetPath, int stride);

void hdf5CreateExtendableInt2DArray(hid_t file_id, const char* datasetPath, int stride);

void hdf5CreateExtendableDouble3DArray(hid_t file_id, const char *datasetPath, int stride1, int stride2);

void hdf5Extend2ndDimensionAndWriteToDouble2DArray(hid_t file_id, const char *datasetPath, const double *buffer);

void hdf5Extend2ndDimensionAndWriteToInt2DArray(hid_t file_id, const char *datasetPath, const int *buffer);

void hdf5CreateOrExtendAndWriteToDouble2DArray(hid_t file_id, const char *parentPath, const char *datasetName, const double *buffer, int stride);

void hdf5CreateOrExtendAndWriteToInt2DArray(hid_t file_id, const char *parentPath, const char *datasetName, const int *buffer, int stride);

void hdf5CreateOrExtendAndWriteToDouble3DArray(hid_t file_id, const char *parentPath, const char *datasetName, const double *buffer, int stride1, int stride2);


#endif
