#include "optimizationOptions.h"
#include "localOptimizationCeres.h"
#include "localOptimizationIpopt.h"
#include "logging.h"
#include "misc.h"
#include <cassert>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <limits>
#include <hdf5.h>
#include <H5Cpp.h>

// Workaround for missing to_string on some systems
namespace patch {
template <typename T> std::string to_string(const T &n) {
    std::ostringstream stm;
    stm << n;
    return stm.str();
}
}



Optimizer *OptimizationOptions::createOptimizer() {
    return optimizerFactory(optimizer);
}

OptimizationOptions *OptimizationOptions::fromHDF5(const char *fileName) {
    hid_t fileId = H5Fopen(fileName, H5F_ACC_RDONLY, H5P_DEFAULT);

    if (fileId < 0) {
        logmessage(
            LOGLVL_CRITICAL,
            "OptimizationOptions::fromHDF5 failed to open HDF5 file '%s'.",
            fileName);
        printBacktrace(20);
    }

    OptimizationOptions *o = fromHDF5(fileId);

    H5Fclose(fileId);

    return o;
}

OptimizationOptions *OptimizationOptions::fromHDF5(hid_t fileId) {
    OptimizationOptions *o = new OptimizationOptions();

    const char *hdf5path = "/optimizationOptions";

    hid_t attributeGroup = H5Gopen1(fileId, hdf5path);
    if(attributeGroup < 0)
        return o;

    //    int numAttributes = H5Aget_num_attrs(attributeGroup);

    H5Aiterate2(attributeGroup, H5_INDEX_NAME, H5_ITER_NATIVE, 0,
                [](hid_t location_id/*in*/,
                const char *attr_name/*in*/,
                const H5A_info_t *ainfo/*in*/,
                void *op_data/*in,out*/) -> herr_t {

        auto *o = static_cast<OptimizationOptions*>(op_data);

        hid_t a = H5Aopen(location_id, attr_name, H5P_DEFAULT);
        hid_t type = H5Aget_type(a);
        hid_t typeClass = H5Tget_class(type);
        hid_t nativeType = H5Tget_native_type(type, H5T_DIR_ASCEND);
        char buf[ainfo->data_size];
        H5Aread(a, nativeType, buf);
        H5Tclose(nativeType);
        H5Aclose(a);

        if (typeClass == H5T_STRING) {
            o->setOption(attr_name, buf);
        } else if (typeClass == H5T_FLOAT) {
            o->setOption(attr_name, *reinterpret_cast<double*>(buf));
        } else if (typeClass == H5T_INTEGER) {
            o->setOption(attr_name, *reinterpret_cast<int*>(buf));
        } else {
            // invalid option type
            abort();
        }

        return 0; // continue
    }, o);


    if (hdf5AttributeExists(fileId, hdf5path, "optimizer")) {
        H5LTget_attribute_int(fileId, hdf5path, "optimizer",
                              (int *)&o->optimizer);
    }

    if (hdf5AttributeExists(fileId, hdf5path, "maxIter")) {
        H5LTget_attribute_int(fileId, hdf5path, "maxIter",
                              &o->maxOptimizerIterations);
    }

    if (hdf5AttributeExists(fileId, hdf5path, "numStarts")) {
        H5LTget_attribute_int(fileId, hdf5path, "numStarts", &o->numStarts);
    }

    if (hdf5AttributeExists(fileId, hdf5path, "retryOptimization")) {
        H5LTget_attribute_int(fileId, hdf5path, "retryOptimization",
                              &o->retryOptimization);
    }

    if (hdf5AttributeExists(fileId, hdf5path, "functionTolerance")) {
        H5LTget_attribute_double(fileId, hdf5path, "functionTolerance",
                                 &o->functionTolerance);
    }

    if (hdf5AttributeExists(fileId, hdf5path, "maxIter")) {
        H5LTget_attribute_int(fileId, hdf5path, "maxIter",
                              &o->maxOptimizerIterations);
    }

    if (hdf5AttributeExists(fileId, hdf5path, "watchdog_shortened_iter_trigger")) {
        H5LTget_attribute_int(fileId, hdf5path, "watchdog_shortened_iter_trigger",
                              &o->watchdog_shortened_iter_trigger);
    }

//    if (hdf5AttributeExists(fileId, hdf5path, "accept_every_trial_step")) {
//        H5LTget_attribute_int(fileId, hdf5path, "accept_every_trial_step",
//                              (int*)&o->accept_every_trial_step);
//    }

    H5Gclose(attributeGroup);
    return o;
}

/**
 * @brief Read pre-generated random starting points from data file.
 * This is used to ensure that the same random starting points
 * are used on different systems where RNG initialization might yield different
 * results.
 * Reads a column from a matrix with random starting points. The size is
 * determined from the dataset.
 * @param fileId
 * @param index Column to read
 * @return The selected starting point or NULL if the dataset did not exist or
 * had less columns than `ìndex`
 */
double *OptimizationOptions::getStartingPoint(hid_t fileId, int index) {
    double *buffer = NULL;

    const char *path = "/optimizationOptions/randomStarts";

    hdf5LockMutex();
    H5_SAVE_ERROR_HANDLER;

    hid_t dataset = H5Dopen2(fileId, path, H5P_DEFAULT);
    if (dataset < 0) {
        logmessage(LOGLVL_DEBUG, "No initial parameters found in %s", path);
        goto freturn;
    }

    {
        // read dimensions
        hid_t dataspace = H5Dget_space(dataset);
        const int ndims = H5Sget_simple_extent_ndims(dataspace);
        assert(ndims == 2);
        hsize_t dims[ndims];
        H5Sget_simple_extent_dims(dataspace, dims, NULL);
        if (dims[1] < (unsigned)index)
            goto freturn;

        logmessage(LOGLVL_INFO, "Reading random initial theta %d from %s",
                   index, path);

        buffer = new double[dims[0]];
        hdf5Read2DDoubleHyperslab(fileId, path, dims[0], 1, 0, index, buffer);
    }

freturn:
    if (H5Eget_num(H5E_DEFAULT)) {
        error("Problem in OptimizationOptions::getStartingPoint\n");
        H5Ewalk2(H5E_DEFAULT, H5E_WALK_DOWNWARD, hdf5ErrorStackWalker_cb, NULL);
    }

    H5_RESTORE_ERROR_HANDLER;
    hdf5UnlockMutex();

    return buffer;
}

std::string OptimizationOptions::toString() {
    std::string s;
    s += "optimizer: " + patch::to_string(optimizer) + "\n";
    s += "maxIter: " + patch::to_string(maxOptimizerIterations) + "\n";
    s += "printToStdout: " + patch::to_string(printToStdout) + "\n";
    s += "numStarts: " + patch::to_string(numStarts) + "\n";
    s += "functionTolerance: " + patch::to_string(functionTolerance) + "\n";
    return s;
}

int OptimizationOptions::getIntOption(std::string key)
{
    return std::stoi(options[key]);
}

double OptimizationOptions::getDoubleOption(std::string key)
{
    return std::stod(options[key]);
}

std::string OptimizationOptions::getStringOption(std::string key)
{
    return options[key];
}

void OptimizationOptions::setOption(std::string key, int value)
{
    options[key] = std::to_string(value);
}

void OptimizationOptions::setOption(std::string key, double value)
{
    std::ostringstream out;
    out << std::setprecision(std::numeric_limits<double>::max_digits10) << value;
    options[key] = out.str();
}

void OptimizationOptions::setOption(std::string key, std::string value)
{
    options[key] = value;
}

Optimizer* optimizerFactory(optimizerEnum optimizer)
{
    switch (optimizer) {
    case OPTIMIZER_IPOPT:
        return new OptimizerIpOpt();
    case OPTIMIZER_CERES:
        return new OptimizerCeres();
    }

    return nullptr;
}
