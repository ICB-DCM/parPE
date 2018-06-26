#include "optimizationResultWriter.h"
#include <parpeVersion.h>

#include <cassert>
#include <cmath>
#include <sstream>
#include <iostream>

#include <hdf5Misc.h>
#include <misc.h>
#include <logging.h>

namespace parpe {

OptimizationResultWriter::OptimizationResultWriter() {}

OptimizationResultWriter::OptimizationResultWriter(hid_t file_id)
{
    auto lock = hdf5MutexGetLock();
    this->file_id = H5Freopen(file_id);

    logParPEVersion();
}


OptimizationResultWriter::OptimizationResultWriter(const std::string &filename,
                                                   bool overwrite) {
    logmessage(LOGLVL_DEBUG, "Writing results to %s.", filename.c_str());
    mkpathConstChar(filename.c_str(), 0755);
    file_id = hdf5CreateFile(filename.c_str(), overwrite);

    if(file_id < 0)
        throw(HDF5Exception());

    logParPEVersion();
}


std::string OptimizationResultWriter::getOptimizationPath() const { return rootPath; }


std::string OptimizationResultWriter::getIterationPath(int iterationIdx) const {
    std::ostringstream ss;
    ss << rootPath << "/" << iterationIdx << "/";

    return ss.str();
}


void OptimizationResultWriter::logParPEVersion() const {
    hdf5EnsureGroupExists(file_id, rootPath.c_str());
    hdf5WriteStringAttribute(file_id, rootPath.c_str(), "PARPE_VERSION",
                             PARPE_VERSION);
}


void OptimizationResultWriter::closeResultHDFFile() {
    auto lock = hdf5MutexGetLock();

    H5_SAVE_ERROR_HANDLER;
    herr_t status = H5Fclose(file_id);

    if (status < 0) {
        error("closeResultHDFFile failed to close HDF5 file.");
        H5Eprint(H5E_DEFAULT, stderr);
    }
    H5_RESTORE_ERROR_HANDLER;
}


void OptimizationResultWriter::logLocalOptimizerObjectiveFunctionEvaluation(gsl::span<const double> parameters, double objectiveFunctionValue,
        gsl::span<const double> objectiveFunctionGradient, int numIterations, int numFunctionCalls,
        double timeElapsedInSeconds) {

    std::string pathStr = getIterationPath(numIterations);
    const char *fullGroupPath = pathStr.c_str();

    hdf5CreateOrExtendAndWriteToDouble2DArray(
                file_id, fullGroupPath, "costFunCost", &objectiveFunctionValue, 1);

    if (!objectiveFunctionGradient.empty()) {
        hdf5CreateOrExtendAndWriteToDouble2DArray(
                    file_id, fullGroupPath, "costFunGradient",
                    objectiveFunctionGradient.data(), objectiveFunctionGradient.size());
    } else if (!parameters.empty()) {
        double dummyGradient[parameters.size()];
        std::fill_n(dummyGradient, parameters.size(), NAN);
        hdf5CreateOrExtendAndWriteToDouble2DArray(file_id, fullGroupPath,
                                                  "costFunGradient",
                                                  dummyGradient, parameters.size());
    }

    if (!parameters.empty())
        hdf5CreateOrExtendAndWriteToDouble2DArray(file_id, fullGroupPath,
                                                  "costFunParameters",
                                                  parameters.data(), parameters.size());

    hdf5CreateOrExtendAndWriteToDouble2DArray(file_id, fullGroupPath,
                                              "costFunWallTimeInSec",
                                              &timeElapsedInSeconds, 1);

    hdf5CreateOrExtendAndWriteToInt2DArray(
                file_id, fullGroupPath, "costFunCallIndex", &numFunctionCalls, 1);

    flushResultWriter();
}


void OptimizationResultWriter::logLocalOptimizerIteration(int numIterations, gsl::span<const double> parameters,
                                                          double objectiveFunctionValue, gsl::span<const double> gradient,
                                                          double timeElapsedInSeconds)
{
    std::string pathStr = getOptimizationPath();
    const char *fullGroupPath = pathStr.c_str();

    hdf5CreateOrExtendAndWriteToDouble2DArray(
                file_id, fullGroupPath, "iterCostFunCost", &objectiveFunctionValue, 1);
    if (!gradient.empty())
        hdf5CreateOrExtendAndWriteToDouble2DArray(file_id, fullGroupPath,
                                                  "iterCostFunGradient",
                                                  gradient.data(), gradient.size());
    if (!parameters.empty())
        hdf5CreateOrExtendAndWriteToDouble2DArray(file_id, fullGroupPath,
                                                  "iterCostFunParameters",
                                                  parameters.data(), parameters.size());
    hdf5CreateOrExtendAndWriteToDouble2DArray(file_id, fullGroupPath,
                                              "iterCostFunWallTimeInSec",
                                              &timeElapsedInSeconds, 1);

    hdf5CreateOrExtendAndWriteToInt2DArray(file_id, fullGroupPath, "iterIndex",
                                           &numIterations, 1);
    /*
    hdf5CreateOrExtendAndWriteToInt2DArray(file_id, fullGroupPath,
                                           "iterIpopt_alg_mod", &alg_mod, 1);
    hdf5CreateOrExtendAndWriteToDouble2DArray(file_id, fullGroupPath,
                                              "iterIpopt_inf_pr", &inf_pr, 1);
    hdf5CreateOrExtendAndWriteToDouble2DArray(file_id, fullGroupPath,
                                              "iterIpopt_inf_du", &inf_du, 1);
    hdf5CreateOrExtendAndWriteToDouble2DArray(file_id, fullGroupPath,
                                              "iterIpopt_mu", &mu, 1);
    hdf5CreateOrExtendAndWriteToDouble2DArray(file_id, fullGroupPath,
                                              "iterIpopt_d_norm", &d_norm, 1);
    hdf5CreateOrExtendAndWriteToDouble2DArray(file_id, fullGroupPath,
                                              "iterIpopt_regularization_size",
                                              &regularization_size, 1);
    hdf5CreateOrExtendAndWriteToDouble2DArray(
        file_id, fullGroupPath, "iterIpopt_alpha_du", &alpha_du, 1);
    hdf5CreateOrExtendAndWriteToDouble2DArray(
        file_id, fullGroupPath, "iterIpopt_alpha_pr", &alpha_pr, 1);
    hdf5CreateOrExtendAndWriteToInt2DArray(
        file_id, fullGroupPath, "iterIpopt_ls_trials", &ls_trials, 1);
*/
    flushResultWriter();
}


void OptimizationResultWriter::starting(gsl::span<double const> initialParameters)
{
    if (!initialParameters.empty()) {
        std::string pathStr = getOptimizationPath();
        const char *fullGroupPath = pathStr.c_str();

        hdf5CreateOrExtendAndWriteToDouble2DArray(file_id, fullGroupPath,
                                                  "initialParameters",
                                                  initialParameters.data(), initialParameters.size());
        flushResultWriter();
    }
}


void OptimizationResultWriter::flushResultWriter() const {
    auto lock = hdf5MutexGetLock();

    H5Fflush(file_id, H5F_SCOPE_LOCAL);
}

void OptimizationResultWriter::saveLocalOptimizerResults(double finalNegLogLikelihood, gsl::span<const double> optimalParameters,
                                                         double masterTime, int exitStatus) const {

    std::string optimPath = getOptimizationPath();
    hdf5EnsureGroupExists(file_id, optimPath.c_str());

    std::string fullGroupPath;
    hsize_t dimensions[1] = {1};

    auto lock = hdf5MutexGetLock();

    fullGroupPath = (optimPath + "/finalCost");
    H5LTmake_dataset(file_id, fullGroupPath.c_str(), 1, dimensions,
                     H5T_NATIVE_DOUBLE, &finalNegLogLikelihood);

    fullGroupPath = (optimPath + "/wallTimeInSec");
    H5LTmake_dataset(file_id, fullGroupPath.c_str(), 1, dimensions,
                     H5T_NATIVE_DOUBLE, &masterTime);

    fullGroupPath = (optimPath + "/exitStatus");
    H5LTmake_dataset(file_id, fullGroupPath.c_str(), 1, dimensions,
                     H5T_NATIVE_INT, &exitStatus);

    if (!optimalParameters.empty()) {
        fullGroupPath = (optimPath + "/finalParameters");
        dimensions[0] = optimalParameters.size();
        H5LTmake_dataset(file_id, fullGroupPath.c_str(), 1, dimensions,
                         H5T_NATIVE_DOUBLE, optimalParameters.data());
    }

    flushResultWriter();
}

void OptimizationResultWriter::setRootPath(const std::string &path)
{
    rootPath = path;
    hdf5EnsureGroupExists(file_id, rootPath.c_str());
}

OptimizationResultWriter::~OptimizationResultWriter() {
    closeResultHDFFile();
}

hid_t OptimizationResultWriter::getFileId() const
{
    return file_id;
}

} // namespace parpe
