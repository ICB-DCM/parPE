#include <parpeoptimization/optimizationResultWriter.h>

#include <parpecommon/parpeVersion.h>
#include <parpecommon/hdf5Misc.h>
#include <parpecommon/misc.h>
#include <parpecommon/logging.h>

#include <cmath>
#include <sstream>
#include <iostream>

#include <hdf5_hl.h>

namespace parpe {

OptimizationResultWriter::OptimizationResultWriter(const H5::H5File file,
                                                   std::string rootPath) :
    rootPath(std::move(rootPath)) {
    auto lock = hdf5MutexGetLock();
    this->file = file;

    hdf5EnsureGroupExists(file, this->rootPath);
}

OptimizationResultWriter::OptimizationResultWriter(const std::string &filename,
                                                   bool overwrite,
                                                   std::string rootPath) :
    rootPath(std::move(rootPath)) {
    logmessage(LOGLVL_DEBUG, "Writing results to %s.", filename.c_str());

    file = hdf5CreateFile(filename.c_str(), overwrite);

    hdf5EnsureGroupExists(file, this->rootPath);
}

OptimizationResultWriter::OptimizationResultWriter(const OptimizationResultWriter &other) :
    rootPath(other.rootPath) {
    auto lock = hdf5MutexGetLock();
    file = other.file;
    hdf5EnsureGroupExists(file, rootPath);
}

const std::string &OptimizationResultWriter::getRootPath() const {
    return rootPath;
}

std::string OptimizationResultWriter::getIterationPath(int iterationIdx) const {
    std::ostringstream ss;
    ss << rootPath << "/iteration/" << iterationIdx << "/";

    return ss.str();
}

void OptimizationResultWriter::logObjectiveFunctionEvaluation(
        gsl::span<const double> parameters,
        double objectiveFunctionValue,
        gsl::span<const double> objectiveFunctionGradient,
        int numIterations,
        int numFunctionCalls,
        double timeElapsedInSeconds)
{

    std::string pathStr = getIterationPath(numIterations);
    const char *fullGroupPath = pathStr.c_str();

    hdf5CreateOrExtendAndWriteToDouble2DArray(
                file.getId(), fullGroupPath, "costFunCost",
                gsl::make_span<const double>(&objectiveFunctionValue, 1));

    if (logGradientEachFunctionEvaluation) {
        if (!objectiveFunctionGradient.empty()) {
            hdf5CreateOrExtendAndWriteToDouble2DArray(
                        file.getId(), fullGroupPath, "costFunGradient",
                        objectiveFunctionGradient);
        } else if (!parameters.empty()) {
            double dummyGradient[parameters.size()];
            std::fill_n(dummyGradient, parameters.size(), NAN);
            hdf5CreateOrExtendAndWriteToDouble2DArray(
                        file.getId(), fullGroupPath, "costFunGradient",
                        gsl::make_span<const double>(dummyGradient,
                                                     parameters.size()));
        }
    }

    if (logParametersEachFunctionEvaluation)
        if (!parameters.empty())
            hdf5CreateOrExtendAndWriteToDouble2DArray(
                        file.getId(), fullGroupPath, "costFunParameters", parameters);

    hdf5CreateOrExtendAndWriteToDouble2DArray(
                file.getId(), fullGroupPath, "costFunWallTimeInSec",
                gsl::make_span<const double>(&timeElapsedInSeconds, 1));

    hdf5CreateOrExtendAndWriteToInt2DArray(
                file.getId(), fullGroupPath, "costFunCallIndex",
                gsl::make_span<const int>(&numFunctionCalls, 1));

    flushResultWriter();
}

void OptimizationResultWriter::logOptimizerIteration(
        int numIterations,
        gsl::span<const double> parameters,
        double objectiveFunctionValue,
        gsl::span<const double> gradient,
        double wallSeconds,
        double cpuSeconds) {
    std::string const& pathStr = getRootPath();
    const char *fullGroupPath = pathStr.c_str();

    hdf5CreateOrExtendAndWriteToDouble2DArray(
                file.getId(), fullGroupPath, "iterCostFunCost",
                gsl::make_span<const double>(&objectiveFunctionValue, 1));

    if (logGradientEachIteration) {
        if (!gradient.empty()) {
            hdf5CreateOrExtendAndWriteToDouble2DArray(
                        file.getId(), fullGroupPath, "iterCostFunGradient",
                        gradient);
        } else if (!parameters.empty()) {
            std::vector<double> nanGradient(parameters.size(), NAN);
            hdf5CreateOrExtendAndWriteToDouble2DArray(
                        file.getId(), fullGroupPath, "iterCostFunGradient",
                        nanGradient);
        }
    }

    if (!parameters.empty()) {
        hdf5CreateOrExtendAndWriteToDouble2DArray(
                    file.getId(), fullGroupPath, "iterCostFunParameters",
                    parameters);
    }

    hdf5CreateOrExtendAndWriteToDouble2DArray(
                file.getId(), fullGroupPath, "iterCostFunWallSec",
                gsl::make_span<const double>(&wallSeconds, 1));

    hdf5CreateOrExtendAndWriteToDouble2DArray(
                file.getId(), fullGroupPath, "iterCostFunCpuSec",
                gsl::make_span(&cpuSeconds, 1));

    hdf5CreateOrExtendAndWriteToInt2DArray(
                file.getId(), fullGroupPath, "iterIndex",
                gsl::make_span<const int>(&numIterations, 1));

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
        std::string const& pathStr = getRootPath();
        const char *fullGroupPath = pathStr.c_str();

        hdf5CreateOrExtendAndWriteToDouble2DArray(
                    file.getId(), fullGroupPath, "initialParameters",
                    initialParameters);
        flushResultWriter();
    }
}

void OptimizationResultWriter::flushResultWriter() const {
    auto lock = hdf5MutexGetLock();

    file.flush(H5F_SCOPE_LOCAL);
}

void OptimizationResultWriter::saveOptimizerResults(
        double finalNegLogLikelihood,
        gsl::span<const double> optimalParameters,
        double wallSec,
        double cpuSec,
        int exitStatus) const {

    std::string const& optimPath = getRootPath();
    hdf5EnsureGroupExists(file, optimPath.c_str());

    std::string fullGroupPath;
    hsize_t dimensions[1] = { 1 };

    auto lock = hdf5MutexGetLock();

    fullGroupPath = (optimPath + "/finalCost");
    H5LTmake_dataset(file.getId(), fullGroupPath.c_str(), 1, dimensions,
                     H5T_NATIVE_DOUBLE, &finalNegLogLikelihood);

    fullGroupPath = (optimPath + "/wallSec");
    H5LTmake_dataset(file.getId(), fullGroupPath.c_str(), 1, dimensions,
                     H5T_NATIVE_DOUBLE, &wallSec);

    fullGroupPath = (optimPath + "/cpuSec");
    H5LTmake_dataset(file.getId(), fullGroupPath.c_str(), 1, dimensions,
                     H5T_NATIVE_DOUBLE, &cpuSec);

    fullGroupPath = (optimPath + "/exitStatus");
    H5LTmake_dataset(file.getId(), fullGroupPath.c_str(), 1, dimensions,
                     H5T_NATIVE_INT, &exitStatus);

    if (!optimalParameters.empty()) {
        fullGroupPath = (optimPath + "/finalParameters");
        dimensions[0] = optimalParameters.size();
        H5LTmake_dataset(file.getId(), fullGroupPath.c_str(), 1, dimensions,
                         H5T_NATIVE_DOUBLE, optimalParameters.data());
    }

    flushResultWriter();
}

void OptimizationResultWriter::setRootPath(const std::string &path) {
    rootPath = path;
    hdf5EnsureGroupExists(file, rootPath.c_str());
}

OptimizationResultWriter::~OptimizationResultWriter() {
}

const H5::H5File &OptimizationResultWriter::getH5File() const {
    return file;
}

} // namespace parpe
