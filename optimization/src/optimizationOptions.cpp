#include "optimizationOptions.h"
#include "localOptimizationCeres.h"
#include "localOptimizationIpopt.h"
#ifdef PARPE_DLIB_ENABLED
#include "localOptimizationDlib.h"
#endif
#ifdef PARPE_TOMS611_ENABLED
#include "localOptimizationToms611.h"
#endif
#ifdef PARPE_FSQP_ENABLED
#include "localOptimizationFsqp.h"
#endif
#include "logging.h"
#include "misc.h"
#include <cassert>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <limits>
#include <hdf5.h>
#include <H5Cpp.h>

namespace parpe {

// Workaround for missing to_string on some systems
namespace patch {
template <typename T> std::string to_string(const T &n) {
    std::ostringstream stm;
    stm << n;
    return stm.str();
}
}


herr_t optimizationOptionsFromAttribute(hid_t location_id/*in*/,
                                        const char *attr_name/*in*/,
                                        const H5A_info_t *ainfo/*in*/,
                                        void *op_data/*in,out*/) {
    // iterate over attributes and add to OptimizationOptions

    auto *o = static_cast<OptimizationOptions*>(op_data);

    hid_t a = H5Aopen(location_id, attr_name, H5P_DEFAULT);
    hid_t type = H5Aget_type(a);
    hid_t typeClass = H5Tget_class(type);
    hid_t nativeType = H5Tget_native_type(type, H5T_DIR_ASCEND);
    char buf[ainfo->data_size + 1]; // +1 for \0
    buf[ainfo->data_size] = '\0';
    H5Aread(a, nativeType, buf);
    H5Tclose(nativeType);
    H5Aclose(a);

    if (typeClass == H5T_STRING) {
        // NOTE: only works for (fixed-length?) ASCII strings, no unicode -> in python use np.string_("bla")
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
}

Optimizer *OptimizationOptions::createOptimizer() const {
    return optimizerFactory(optimizer);
}

std::unique_ptr<OptimizationOptions> OptimizationOptions::fromHDF5(const char *fileName) {
    H5::H5File file;
    try {
        file = H5::H5File(fileName, H5F_ACC_RDONLY);
    } catch (...) {
        throw HDF5Exception(
                    "OptimizationOptions::fromHDF5 failed to open HDF5 file '%s'.",
                    fileName);
    }

    return fromHDF5(file.getId());
}

std::unique_ptr<OptimizationOptions> OptimizationOptions::fromHDF5(hid_t fileId) {
    auto o = std::make_unique<OptimizationOptions>();

    const char *hdf5path = "/optimizationOptions";


    if (hdf5AttributeExists(fileId, hdf5path, "optimizer")) {
        H5LTget_attribute_int(fileId, hdf5path, "optimizer",
                              (int *)&o->optimizer);
    }

    if (hdf5AttributeExists(fileId, hdf5path, "numStarts")) {
        H5LTget_attribute_int(fileId, hdf5path, "numStarts", &o->numStarts);
    }

    if (hdf5AttributeExists(fileId, hdf5path, "retryOptimization")) {
        H5LTget_attribute_int(fileId, hdf5path, "retryOptimization",
                              &o->retryOptimization);
    }

    if (hdf5AttributeExists(fileId, hdf5path, "multistartsInParallel")) {
        H5LTget_attribute_int(fileId, hdf5path, "multistartsInParallel",
                              &o->multistartsInParallel);
    }

    if (hdf5AttributeExists(fileId, hdf5path, "maxIter")) {
        // this value is overwritten by any optimizer-specific configuration
        H5LTget_attribute_int(fileId, hdf5path, "maxIter", &o->maxOptimizerIterations);
    }

    std::string optimizerPath;

    switch(o->optimizer) {
    case OPTIMIZER_CERES:
        optimizerPath = std::string(hdf5path) + "/ceres";
        break;
    case OPTIMIZER_TOMS611:
        optimizerPath = std::string(hdf5path) + "/toms611";
        break;
    case OPTIMIZER_IPOPT:
    default:
        optimizerPath = std::string(hdf5path) + "/ipopt";
    }

    if(hdf5GroupExists(fileId, optimizerPath.c_str())) {
        hid_t attributeGroup = H5Gopen1(fileId, optimizerPath.c_str());
        if(attributeGroup < 0)
            return o;

        H5Aiterate2(attributeGroup, H5_INDEX_NAME, H5_ITER_NATIVE, 0,
                    optimizationOptionsFromAttribute, o.get());

        H5Gclose(attributeGroup);
    }
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
std::vector<double> OptimizationOptions::getStartingPoint(hid_t fileId, int index) {
    std::vector<double> startingPoint;

    const char *path = "/optimizationOptions/randomStarts";

    hdf5LockMutex();
    H5_SAVE_ERROR_HANDLER;

    hid_t dataset;
    if (!hdf5DatasetExists(fileId, path) || (dataset = H5Dopen2(fileId, path, H5P_DEFAULT)) < 0) {
        logmessage(LOGLVL_DEBUG, "No initial parameters found in %s", path);
        H5Eclear1();
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

        startingPoint.resize(dims[0]);
        hdf5Read2DDoubleHyperslab(fileId, path, dims[0], 1, 0, index, startingPoint.data());
    }

freturn:
    if (H5Eget_num(H5E_DEFAULT)) {
        error("Problem in OptimizationOptions::getStartingPoint\n");
        H5Ewalk2(H5E_DEFAULT, H5E_WALK_DOWNWARD, hdf5ErrorStackWalker_cb, NULL);
    }

    H5_RESTORE_ERROR_HANDLER;
    hdf5UnlockMutex();

    return startingPoint;
}

std::string OptimizationOptions::toString() {
    std::string s;
    s += "optimizer: " + patch::to_string(optimizer) + "\n";
    s += "maxIter: " + patch::to_string(maxOptimizerIterations) + "\n";
    s += "printToStdout: " + patch::to_string(printToStdout) + "\n";
    s += "numStarts: " + patch::to_string(numStarts) + "\n";
    s += "\n";

    for_each<std::string&>([](const std::pair<const std::string, const std::string> pair, std::string &out){
        out = out + pair.first + ": " + pair.second + "\n";
    }, s);

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
    case OPTIMIZER_DLIB:
#ifdef PARPE_DLIB_ENABLED
        return new OptimizerDlibLineSearch();
#else
        return nullptr;
#endif
    case OPTIMIZER_TOMS611:
#ifdef PARPE_TOMS611_ENABLED
        return new OptimizerToms611TrustRegionSumsl();
#else
        return nullptr;
#endif
    case OPTIMIZER_FSQP:
#ifdef PARPE_FSQP_ENABLED
        return new OptimizerFsqp();
#else
        return nullptr;
#endif

    }

    return nullptr;
}

} // namespace parpe
