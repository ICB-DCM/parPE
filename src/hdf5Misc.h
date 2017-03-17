#ifndef HDF5_MISC_H
#define HDF5_MISC_H

#include <hdf5.h>
#include <hdf5_hl.h>

#include <pthread.h>

// mutex for **ALL** HDF5 library calls; read and write; any file(?)
pthread_mutex_t mutexHDF;

void initHDF5Mutex();


#define H5_SAVE_ERROR_HANDLER   herr_t (*old_func)(void*); \
                                void *old_client_data; \
                                H5Eget_auto1(&old_func, &old_client_data); \
                                H5Eset_auto1(NULL, NULL);

#define H5_RESTORE_ERROR_HANDLER H5Eset_auto1(old_func, old_client_data);

herr_t hdf5ErrorStackWalker_cb(unsigned int n, const H5E_error_t *err_desc, void *client_data); // TODO: also use for resultwriter


// malloc version of ami_hdf5.cpp
void getDoubleArrayAttributeC(hid_t file_id, const char* optionsObject, const char* attributeName, double **destination, hsize_t *length);
void getIntArrayAttributeC(hid_t file_id, const char* optionsObject, const char* attributeName, int **destination, hsize_t *length);

#endif
