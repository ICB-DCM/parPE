#include "optimizationOptions.h"
#include "logging.h"
#include "misc.h"
#include <string>

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

    std::cout<<o->toString();

    return o;
}

std::string OptimizationOptions::toString()
{
    std::string s;
    s += "optimizer: " + std::to_string(optimizer) + "\n";
    s += "maxIter: " + std::to_string(maxOptimizerIterations) + "\n";
    s += "printToStdout: " + std::to_string(printToStdout) + "\n";
    s += "numStarts: " + std::to_string(numStarts) + "\n";
    s += "functionTolerance: " + std::to_string(functionTolerance) + "\n";
    return s;
}
