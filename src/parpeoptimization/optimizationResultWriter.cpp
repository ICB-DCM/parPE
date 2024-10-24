#include <parpeoptimization/optimizationResultWriter.h>

#include <parpecommon/hdf5Misc.h>
#include <parpecommon/logging.h>
#include <parpecommon/misc.h>
#include <parpecommon/parpeVersion.h>

#include <cmath>
#include <iostream>
#include <sstream>

#include <hdf5_hl.h>

namespace parpe {

OptimizationResultWriter::OptimizationResultWriter(
    const H5::H5File& file,
    std::string rootPath)
    : rootPath(std::move(rootPath)) {
    [[maybe_unused]] auto lock = hdf5MutexGetLock();
    this->file = file;

    hdf5EnsureGroupExists(file, this->rootPath);
}

OptimizationResultWriter::OptimizationResultWriter(
    std::string const& filename,
    bool overwrite,
    std::string rootPath)
    : rootPath(std::move(rootPath)) {
    logmessage(loglevel::debug, "Writing results to %s.", filename.c_str());

    file = hdf5CreateFile(filename, overwrite);

    hdf5EnsureGroupExists(file, this->rootPath);
}

OptimizationResultWriter::OptimizationResultWriter(
    OptimizationResultWriter const& other)
    : rootPath(other.rootPath) {
    [[maybe_unused]] auto lock = hdf5MutexGetLock();
    file = other.file;
    hdf5EnsureGroupExists(file, rootPath);
}

OptimizationResultWriter::~OptimizationResultWriter() {
    [[maybe_unused]] auto lock = hdf5MutexGetLock();
    file.close();
}

std::string const& OptimizationResultWriter::getRootPath() const {
    return rootPath;
}

std::string OptimizationResultWriter::getIterationPath(int iterationIdx) const {
    std::ostringstream ss;
    ss << rootPath << "/iteration/" << iterationIdx << "/";

    return ss.str();
}

void OptimizationResultWriter::logObjectiveFunctionEvaluation(
    gsl::span<double const> parameters,
    double objectiveFunctionValue,
    gsl::span<double const> objectiveFunctionGradient,
    int numIterations,
    int numFunctionCalls,
    double timeElapsedInSeconds) {
    [[maybe_unused]] auto lock = hdf5MutexGetLock();

    std::string pathStr = getIterationPath(numIterations);
    char const* fullGroupPath = pathStr.c_str();

    hdf5CreateOrExtendAndWriteToDouble2DArray(
        file,
        fullGroupPath,
        "costFunCost",
        gsl::make_span<double const>(&objectiveFunctionValue, 1));

    if (logGradientEachFunctionEvaluation) {
        if (!objectiveFunctionGradient.empty()) {
            hdf5CreateOrExtendAndWriteToDouble2DArray(
                file,
                fullGroupPath,
                "costFunGradient",
                objectiveFunctionGradient);
        } else if (!parameters.empty()) {
            double dummyGradient[parameters.size()];
            std::fill_n(dummyGradient, parameters.size(), NAN);
            hdf5CreateOrExtendAndWriteToDouble2DArray(
                file,
                fullGroupPath,
                "costFunGradient",
                gsl::make_span<double const>(dummyGradient, parameters.size()));
        }
    }

    if (logParametersEachFunctionEvaluation)
        if (!parameters.empty())
            hdf5CreateOrExtendAndWriteToDouble2DArray(
                file, fullGroupPath, "costFunParameters", parameters);

    hdf5CreateOrExtendAndWriteToDouble2DArray(
        file,
        fullGroupPath,
        "costFunWallTimeInSec",
        gsl::make_span<double const>(&timeElapsedInSeconds, 1));

    hdf5CreateOrExtendAndWriteToInt2DArray(
        file,
        fullGroupPath,
        "costFunCallIndex",
        gsl::make_span<int const>(&numFunctionCalls, 1));

    flushResultWriter();
}

void OptimizationResultWriter::logOptimizerIteration(
    int numIterations,
    gsl::span<double const> parameters,
    double objectiveFunctionValue,
    gsl::span<double const> gradient,
    double wallSeconds,
    double cpuSeconds) {
    std::string const& pathStr = getRootPath();
    char const* fullGroupPath = pathStr.c_str();

    [[maybe_unused]] auto lock = hdf5MutexGetLock();

    hdf5CreateOrExtendAndWriteToDouble2DArray(
        file,
        fullGroupPath,
        "iterCostFunCost",
        gsl::make_span<double const>(&objectiveFunctionValue, 1));

    if (logGradientEachIteration) {
        if (!gradient.empty()) {
            hdf5CreateOrExtendAndWriteToDouble2DArray(
                file, fullGroupPath, "iterCostFunGradient", gradient);
        } else if (!parameters.empty()) {
            std::vector<double> nanGradient(parameters.size(), NAN);
            hdf5CreateOrExtendAndWriteToDouble2DArray(
                file, fullGroupPath, "iterCostFunGradient", nanGradient);
        }
    }

    if (!parameters.empty()) {
        hdf5CreateOrExtendAndWriteToDouble2DArray(
            file, fullGroupPath, "iterCostFunParameters", parameters);
    }

    hdf5CreateOrExtendAndWriteToDouble2DArray(
        file,
        fullGroupPath,
        "iterCostFunWallSec",
        gsl::make_span<double const>(&wallSeconds, 1));

    hdf5CreateOrExtendAndWriteToDouble2DArray(
        file,
        fullGroupPath,
        "iterCostFunCpuSec",
        gsl::make_span(&cpuSeconds, 1));

    hdf5CreateOrExtendAndWriteToInt2DArray(
        file,
        fullGroupPath,
        "iterIndex",
        gsl::make_span<int const>(&numIterations, 1));

    flushResultWriter();
}

void OptimizationResultWriter::setLoggingEachIteration(bool logGradient) {
    logGradientEachIteration = logGradient;
}

void OptimizationResultWriter::setLoggingEachFunctionEvaluation(
    bool logGradient,
    bool logParameters) {
    logGradientEachFunctionEvaluation = logGradient;
    logParametersEachFunctionEvaluation = logParameters;
}

void OptimizationResultWriter::starting(
    gsl::span<double const> initialParameters) {
    if (!initialParameters.empty()) {
        auto const& root_path = getRootPath();

        hdf5CreateOrExtendAndWriteToDouble2DArray(
            file, root_path, "initialParameters", initialParameters);
        flushResultWriter();
    }
}

void OptimizationResultWriter::flushResultWriter() const {
    [[maybe_unused]] auto lock = hdf5MutexGetLock();

    file.flush(H5F_SCOPE_LOCAL);
}

void OptimizationResultWriter::saveOptimizerResults(
    double finalNegLogLikelihood,
    gsl::span<double const> optimalParameters,
    double wallSec,
    double cpuSec,
    int exitStatus) const {

    std::string const& optimPath = getRootPath();
    hdf5EnsureGroupExists(file, optimPath);

    std::string fullGroupPath;
    hsize_t dimensions[1] = {1};

    [[maybe_unused]] auto lock = hdf5MutexGetLock();

    fullGroupPath = (optimPath + "/finalCost");
    H5LTmake_dataset(
        file.getId(),
        fullGroupPath.c_str(),
        1,
        dimensions,
        H5T_NATIVE_DOUBLE,
        &finalNegLogLikelihood);

    fullGroupPath = (optimPath + "/wallSec");
    H5LTmake_dataset(
        file.getId(),
        fullGroupPath.c_str(),
        1,
        dimensions,
        H5T_NATIVE_DOUBLE,
        &wallSec);

    fullGroupPath = (optimPath + "/cpuSec");
    H5LTmake_dataset(
        file.getId(),
        fullGroupPath.c_str(),
        1,
        dimensions,
        H5T_NATIVE_DOUBLE,
        &cpuSec);

    fullGroupPath = (optimPath + "/exitStatus");
    H5LTmake_dataset(
        file.getId(),
        fullGroupPath.c_str(),
        1,
        dimensions,
        H5T_NATIVE_INT,
        &exitStatus);

    if (!optimalParameters.empty()) {
        fullGroupPath = (optimPath + "/finalParameters");
        dimensions[0] = optimalParameters.size();
        H5LTmake_dataset(
            file.getId(),
            fullGroupPath.c_str(),
            1,
            dimensions,
            H5T_NATIVE_DOUBLE,
            optimalParameters.data());
    }

    flushResultWriter();
}

void OptimizationResultWriter::setRootPath(std::string const& path) {
    rootPath = path;
    hdf5EnsureGroupExists(file, rootPath);
}

const H5::H5File& OptimizationResultWriter::getH5File() const { return file; }

} // namespace parpe
