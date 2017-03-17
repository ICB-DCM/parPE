#include "hdf5Misc.h"

#include "logging.h"
#include "misc.h"
#include <assert.h>
#include <unistd.h>
#include <stdlib.h>

// mutex for **ALL** HDF5 library calls; read and write; any file(?)
static pthread_mutex_t mutexHDF;

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

void getDoubleArrayAttributeC(hid_t file_id, const char* optionsObject, const char* attributeName, double **destination, hsize_t *length) {
    H5T_class_t type_class;
    size_t type_size;
    herr_t status;

    pthread_mutex_lock(&mutexHDF);

    status = H5LTget_attribute_info(file_id, optionsObject, attributeName, length, &type_class, &type_size);
    if(status < 0) {
        fprintf(stderr, "Error in getDoubleArrayAttributeC: Cannot read attribute '%s' of '%s'\n", attributeName, optionsObject);
        printBacktrace(10);
    }

    *destination = malloc(sizeof(double) * (*length)); // vs. type_size
    status = H5LTget_attribute_double(file_id, optionsObject, attributeName, *destination);
    if(status < 0)
        fprintf(stderr, "Error in getDoubleArrayAttributeC: Cannot read attribute '%s' of '%s'\n", attributeName, optionsObject);

    pthread_mutex_unlock(&mutexHDF);
}

void getIntArrayAttributeC(hid_t file_id, const char* optionsObject, const char* attributeName, int **destination, hsize_t *length) {
    H5T_class_t type_class;
    size_t type_size;

    pthread_mutex_lock(&mutexHDF);

    H5LTget_attribute_info(file_id, optionsObject, attributeName, length, &type_class, &type_size);
#ifdef AMI_HDF5_H_DEBUG
    printf("%s: %d: ", attributeName, *length);
#endif
    *destination = (int*) malloc(sizeof(int) * (*length));
    H5LTget_attribute_int(file_id, optionsObject, attributeName, *destination);

    pthread_mutex_unlock(&mutexHDF);
}

