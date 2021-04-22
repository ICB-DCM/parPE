#include <parpecommon/parpeConfig.h>
#include <parpeoptimization/optimizationOptions.h>

#ifdef PARPE_ENABLE_CERES
#include <parpeoptimization/localOptimizationCeres.h>
#endif

#ifdef PARPE_ENABLE_IPOPT
#include <parpeoptimization/localOptimizationIpopt.h>
#endif

#ifdef PARPE_ENABLE_DLIB
#include <parpeoptimization/localOptimizationDlib.h>
#endif

#ifdef PARPE_ENABLE_TOMS611
#include <parpeoptimization/localOptimizationToms611.h>
#endif

#ifdef PARPE_ENABLE_FSQP
#include <parpeoptimization/localOptimizationFsqp.h>
#endif

#include <parpecommon/logging.h>
#include <parpecommon/misc.h>
#include <parpecommon/parpeException.h>

#include <cassert>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <limits>
#include <utility>

#include <H5Cpp.h>

namespace parpe {

// Workaround for missing to_string on some systems
namespace patch {
template <typename T> std::string to_string(const T &n) {
    std::ostringstream stm;
    stm << n;
    return stm.str();
}
} // namespace patch


void optimizationOptionsFromAttribute(H5::H5Object& loc,
                                      const H5std_string attr_name,
                                      void *op_data) {
    // iterate over attributes and add to OptimizationOptions

    auto *o = static_cast<OptimizationOptions*>(op_data);

    auto a = loc.openAttribute(attr_name);
    auto type = a.getDataType();
    auto typeClass = a.getTypeClass();
    H5A_info_t ainfo;
    H5Aget_info(a.getId(), &ainfo);
    char buf[ainfo.data_size + 1]; // +1 for \0
    buf[ainfo.data_size] = '\0';
    auto nativeType = H5Tget_native_type(type.getId(), H5T_DIR_ASCEND);
    a.read(nativeType, buf);
    H5Tclose(nativeType);


    if (typeClass == H5T_STRING) {
        // NOTE: only works for (fixed-length?) ASCII strings, no unicode
        // -> in python use np.string_("bla")
        o->setOption(attr_name, buf);
    } else if (typeClass == H5T_FLOAT) {
        o->setOption(attr_name, *reinterpret_cast<double*>(buf));
    } else if (typeClass == H5T_INTEGER) {
        o->setOption(attr_name, *reinterpret_cast<int*>(buf));
    } else {
        // invalid option type
        abort();
    }
}

std::unique_ptr<Optimizer> OptimizationOptions::createOptimizer() const {
    return optimizerFactory(optimizer);
}

std::unique_ptr<OptimizationOptions> OptimizationOptions::fromHDF5(const std::string &fileName) {
    return fromHDF5(hdf5OpenForReading(fileName));
}

std::unique_ptr<OptimizationOptions> OptimizationOptions::fromHDF5(const H5::H5File &file, std::string const& path) {
    auto o = std::make_unique<OptimizationOptions>();

    const char *hdf5path = path.c_str();
    auto fileId = file.getId();
    [[maybe_unused]] auto lock = hdf5MutexGetLock();

    if (hdf5AttributeExists(file, path, "optimizer")) {
        int buffer;
        H5LTget_attribute_int(fileId, hdf5path, "optimizer", &buffer);
        o->optimizer = static_cast<parpe::optimizerName>(buffer);
    }

    if (hdf5AttributeExists(file, hdf5path, "numStarts")) {
        H5LTget_attribute_int(fileId, hdf5path, "numStarts", &o->numStarts);
    }

    if (hdf5AttributeExists(file, hdf5path, "retryOptimization")) {
        H5LTget_attribute_int(fileId, hdf5path, "retryOptimization",
                              &o->retryOptimization);
    }

    if (hdf5AttributeExists(file, hdf5path, "hierarchicalOptimization")) {
        H5LTget_attribute_int(fileId, hdf5path, "hierarchicalOptimization",
                              &o->hierarchicalOptimization);
    }

    if (hdf5AttributeExists(file, hdf5path, "multistartsInParallel")) {
        H5LTget_attribute_int(fileId, hdf5path, "multistartsInParallel",
                              &o->multistartsInParallel);
    }

    if (hdf5AttributeExists(file, hdf5path, "maxIter")) {
        // this value is overwritten by any optimizer-specific configuration
        H5LTget_attribute_int(fileId, hdf5path, "maxIter", &o->maxOptimizerIterations);
    }

    std::string optimizerPath;

    switch(o->optimizer) {
    case optimizerName::OPTIMIZER_CERES:
        optimizerPath = std::string(hdf5path) + "/ceres";
        break;
    case optimizerName::OPTIMIZER_TOMS611:
        optimizerPath = std::string(hdf5path) + "/toms611";
        break;
    case optimizerName::OPTIMIZER_MINIBATCH_1:
        optimizerPath = std::string(hdf5path) + "/minibatch";
        break;
    case optimizerName::OPTIMIZER_IPOPT:
    default:
        optimizerPath = std::string(hdf5path) + "/ipopt";
    }

    if(hdf5GroupExists(file, optimizerPath)) {
        auto group = file.openGroup(optimizerPath);
        group.iterateAttrs(optimizationOptionsFromAttribute, nullptr, o.get());
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
std::vector<double> OptimizationOptions::getStartingPoint(H5::H5File const& file,
                                                          int index) {
    std::vector<double> startingPoint;

    const char *path = "/optimizationOptions/randomStarts";

    [[maybe_unused]] auto lock = hdf5MutexGetLock();

    if (!file.nameExists(path)) {
        logmessage(LOGLVL_DEBUG, "No initial parameters found in %s", path);
        return startingPoint;
    }

    H5_SAVE_ERROR_HANDLER;
    try {
        auto dataset = file.openDataSet(path);
        // read dimensions
        auto dataspace = dataset.getSpace();
        const int ndims = dataspace.getSimpleExtentNdims();
        Expects(ndims == 2);
        hsize_t dims[ndims];
        dataspace.getSimpleExtentDims(dims);
        if (dims[1] < static_cast<hsize_t>(index)) {
            logmessage(LOGLVL_ERROR,
                       "Requested starting point index %d out of bounds (%d)",
                       index, static_cast<int>(dims[1]));
            return startingPoint;
        }

        logmessage(LOGLVL_INFO, "Reading random initial theta %d from %s",
                   index, path);

        startingPoint.resize(dims[0]);
        hdf5Read2DDoubleHyperslab(file, path, dims[0], 1, 0, index,
                startingPoint);


    }  catch (H5::Exception const&) {
        if (H5Eget_num(H5E_DEFAULT)) {
            error("Problem in OptimizationOptions::getStartingPoint\n");
            H5Ewalk2(H5E_DEFAULT, H5E_WALK_DOWNWARD, hdf5ErrorStackWalker_cb, nullptr);
        }
    }
    H5_RESTORE_ERROR_HANDLER;

    return startingPoint;
}

std::string OptimizationOptions::toString() {
    std::string s;
    s += "optimizer: " + patch::to_string(static_cast<int>(optimizer)) + "\n";
    s += "maxIter: " + patch::to_string(maxOptimizerIterations) + "\n";
    s += "printToStdout: " + patch::to_string(printToStdout) + "\n";
    s += "numStarts: " + patch::to_string(numStarts) + "\n";
    s += "\n";

    for_each<std::string&>(
                [](const std::pair<const std::string, const std::string> pair,
                std::string &out)
    {
        out = out + pair.first + ": " + pair.second + "\n";
    }, s);

    return s;
}

int OptimizationOptions::getIntOption(std::string const& key)
{
    return std::stoi(options[key]);
}

double OptimizationOptions::getDoubleOption(std::string const& key)
{
    return std::stod(options[key]);
}

std::string OptimizationOptions::getStringOption(const std::string& key)
{
    return options[key];
}

void OptimizationOptions::setOption(std::string const& key, int value)
{
    options[key] = std::to_string(value);
}

void OptimizationOptions::setOption(std::string const& key, double value)
{
    std::ostringstream out;
    out << std::setprecision(std::numeric_limits<double>::max_digits10) << value;
    options[key] = out.str();
}

void OptimizationOptions::setOption(const std::string& key, std::string value)
{
    options[key] = std::move(value);
}

std::unique_ptr<Optimizer> optimizerFactory(optimizerName optimizer)
{
    switch (optimizer) {
    case optimizerName::OPTIMIZER_IPOPT:
#ifdef PARPE_ENABLE_IPOPT
        return std::make_unique<OptimizerIpOpt>();
#else
        return nullptr;
#endif
    case optimizerName::OPTIMIZER_CERES:
#ifdef PARPE_ENABLE_CERES
        return std::make_unique<OptimizerCeres>();
#else
        return nullptr;
#endif
    case optimizerName::OPTIMIZER_DLIB:
#ifdef PARPE_ENABLE_DLIB
        return std::make_unique<OptimizerDlibLineSearch>();
#else
        return nullptr;
#endif
    case optimizerName::OPTIMIZER_TOMS611:
#ifdef PARPE_ENABLE_TOMS611
        return std::make_unique<OptimizerToms611TrustRegionSumsl>();
#else
        return nullptr;
#endif
    case optimizerName::OPTIMIZER_FSQP:
#ifdef PARPE_ENABLE_FSQP
        return std::make_unique<OptimizerFsqp>();
#else
        return nullptr;
#endif
    case optimizerName::OPTIMIZER_MINIBATCH_1:
        throw ParPEException("optimizerFactory() cannot be used with "
                             "mini-batch optimizer.");
    }

    return nullptr;
}


void printAvailableOptimizers(std::string const& prefix)
{
    optimizerName optimizer {optimizerName::OPTIMIZER_IPOPT};

    // Note: Keep fall-through switch statement, so compiler will warn us about
    // any addition to optimizerName
    switch (optimizer) {
    case optimizerName::OPTIMIZER_IPOPT:
#ifdef PARPE_ENABLE_IPOPT
        std::cout<<prefix<<std::left<<std::setw(22)<<"OPTIMIZER_IPOPT\t"
                <<static_cast<int>(optimizerName::OPTIMIZER_IPOPT)
               <<" enabled\n";
#else
        std::cout<<prefix<<std::left<<std::setw(22)<<"OPTIMIZER_IPOPT"
                <<static_cast<int>(optimizerName::OPTIMIZER_IPOPT)
               <<" disabled\n";
#endif
        [[fallthrough]];
    case optimizerName::OPTIMIZER_CERES:
#ifdef PARPE_ENABLE_CERES
        std::cout<<prefix<<std::left<<std::setw(22)<<"OPTIMIZER_CERES"
                <<static_cast<int>(optimizerName::OPTIMIZER_CERES)
               <<" enabled\n";
#else
        std::cout<<prefix<<std::left<<std::setw(22)<<"OPTIMIZER_CERES"
                <<static_cast<int>(optimizerName::OPTIMIZER_CERES)
               <<" disabled\n";
#endif
        [[fallthrough]];
    case optimizerName::OPTIMIZER_DLIB:
#ifdef PARPE_ENABLE_DLIB
        std::cout<<prefix<<std::left<<std::setw(22)<<"OPTIMIZER_DLIB"
                <<static_cast<int>(optimizerName::OPTIMIZER_DLIB)
               <<" enabled\n";
#else
        std::cout<<prefix<<std::left<<std::setw(22)<<"OPTIMIZER_DLIB"
                <<static_cast<int>(optimizerName::OPTIMIZER_DLIB)
               <<" disabled\n";
#endif
        [[fallthrough]];
    case optimizerName::OPTIMIZER_TOMS611:
#ifdef PARPE_ENABLE_TOMS611
        std::cout<<prefix<<std::left<<std::setw(22)<<"OPTIMIZER_TOMS611"
                <<static_cast<int>(optimizerName::OPTIMIZER_TOMS611)
               <<" enabled\n";
#else
        std::cout<<prefix<<std::left<<std::setw(22)<<"OPTIMIZER_TOMS611"
                <<static_cast<int>(optimizerName::OPTIMIZER_TOMS611)
               <<" disabled\n";
#endif
        [[fallthrough]];
    case optimizerName::OPTIMIZER_FSQP:
#ifdef PARPE_ENABLE_FSQP
        std::cout<<prefix<<std::left<<std::setw(22)<<"OPTIMIZER_FSQP"
                <<static_cast<int>(optimizerName::OPTIMIZER_FSQP)
               <<" enabled\n";
#else
        std::cout<<prefix<<std::left<<std::setw(22)<<"OPTIMIZER_FSQP"
                <<static_cast<int>(optimizerName::OPTIMIZER_FSQP)
               <<" disabled\n";
#endif
        [[fallthrough]];
    case optimizerName::OPTIMIZER_MINIBATCH_1:
        std::cout<<prefix<<std::left<<std::setw(22)<<"OPTIMIZER_MINIBATCH_1"
                <<static_cast<int>(optimizerName::OPTIMIZER_MINIBATCH_1)
               <<" enabled\n";
    }
}


} // namespace parpe
