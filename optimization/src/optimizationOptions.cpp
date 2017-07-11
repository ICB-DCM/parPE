#include "optimizationOptions.h"
#include "logging.h"
#include "misc.h"
#include <iostream>
#include <sstream>
#include <cassert>

// Workaround for missing to_string on some systems
namespace patch
{
    template < typename T > std::string to_string( const T& n )
    {
        std::ostringstream stm ;
        stm << n ;
        return stm.str() ;
    }
}

OptimizationOptions::OptimizationOptions()
{
    optimizer = OPTIMIZER_IPOPT;
    logFile = NULL;
    printToStdout = true;
    maxOptimizerIterations = 1000;
    numStarts = 1;
    functionTolerance = 1e-18;
}

OptimizationOptions *OptimizationOptions::fromHDF5(const char *fileName)
{
    hid_t fileId = H5Fopen(fileName, H5F_ACC_RDONLY, H5P_DEFAULT);

    if(fileId < 0) {
        logmessage(LOGLVL_CRITICAL, "OptimizationOptions::fromHDF5 failed to open HDF5 file '%s'.", fileName);
        printBacktrace(20);
    }

    OptimizationOptions *o = fromHDF5(fileId);

    H5Fclose(fileId);

    return o;
}

OptimizationOptions *OptimizationOptions::fromHDF5(hid_t fileId)
{
    const char *hdf5path = "optimizationOptions";

    OptimizationOptions *o = new OptimizationOptions();

    if(hdf5AttributeExists(fileId, hdf5path, "optimizer")) {
        H5LTget_attribute_int(fileId, hdf5path, "optimizer", (int *) &o->optimizer);
    }

    if(hdf5AttributeExists(fileId, hdf5path, "maxIter")) {
        H5LTget_attribute_int(fileId, hdf5path, "maxIter", &o->maxOptimizerIterations);
    }

    if(hdf5AttributeExists(fileId, hdf5path, "numStarts")) {
        H5LTget_attribute_int(fileId, hdf5path, "numStarts", &o->numStarts);
    }

    if(hdf5AttributeExists(fileId, hdf5path, "retryOptimization")) {
        H5LTget_attribute_int(fileId, hdf5path, "retryOptimization", &o->retryOptimization);
    }

    if(hdf5AttributeExists(fileId, hdf5path, "functionTolerance")) {
        H5LTget_attribute_double(fileId, hdf5path, "functionTolerance", &o->functionTolerance);
    }

    return o;
}

/**
 * @brief Reads a column from a matrix with random starting points. The size is determined from the dataset.
 * @param fileId
 * @param index Column to read
 * @return The selected starting point or NULL if the dataset did not exist or had less columns than `ìndex`
 */
double *OptimizationOptions::getStartingPoint(hid_t fileId, int index)
{
    const char *path = "/randomstarts";
    logmessage(LOGLVL_INFO, "Reading random initial theta %d from %s", index, path);

    hdf5LockMutex();
    H5_SAVE_ERROR_HANDLER;

    hid_t dataset   = H5Dopen2(fileId, path, H5P_DEFAULT);
    if(dataset < 0)
        return NULL;

    hid_t dataspace = H5Dget_space(dataset);
    const int ndims = H5Sget_simple_extent_ndims(dataspace);
    assert(ndims == 2);
    hsize_t dims[ndims];
    H5Sget_simple_extent_dims(dataspace, dims, NULL);
    if(dims[1] < index)
        return NULL;

    double *buffer = new double[dims[0]];
    hdf5Read2DDoubleHyperslab(fileId, path, dims[0], 1, 0, index, buffer);

    if(H5Eget_num(H5E_DEFAULT)) {
        error("Problem in getRandomInitialThetaFromFile\n");
        H5Ewalk2(H5E_DEFAULT, H5E_WALK_DOWNWARD, hdf5ErrorStackWalker_cb, NULL);
    }

    H5_RESTORE_ERROR_HANDLER;

    hdf5UnlockMutex();

    return buffer;
}

std::string OptimizationOptions::toString()
{
    std::string s;
    s += "optimizer: " + patch::to_string(optimizer) + "\n";
    s += "maxIter: " + patch::to_string(maxOptimizerIterations) + "\n";
    s += "printToStdout: " + patch::to_string(printToStdout) + "\n";
    s += "numStarts: " + patch::to_string(numStarts) + "\n";
    s += "functionTolerance: " + patch::to_string(functionTolerance) + "\n";
    return s;
}
