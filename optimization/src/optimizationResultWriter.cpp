#include "optimizationResultWriter.h"
#include "hdf5Misc.h"
#include "misc.h"
#include <cassert>
#include <cmath>
#include <logging.h>
#include <sstream>
#include <iostream>

namespace parpe {

OptimizationResultWriter::OptimizationResultWriter() {}

OptimizationResultWriter::OptimizationResultWriter(hid_t file_id)
    : file_id(H5Freopen(file_id)) {
    logParPEVersion();
}


OptimizationResultWriter::OptimizationResultWriter(const std::string &filename,
                                                   bool overwrite) {
    logmessage(LOGLVL_DEBUG, "Writing results to %s.", filename.c_str());
    mkpathConstChar(filename.c_str(), 0755);
    file_id = hdf5CreateFile(filename.c_str(), overwrite);

    if(file_id < 0)
        throw(HDF5Exception());

    hdf5EnsureGroupExists(file_id, rootPath.c_str());
    logParPEVersion();
}


std::string OptimizationResultWriter::getOptimizationPath() { return rootPath; }


std::string OptimizationResultWriter::getIterationPath(int iterationIdx) {
    std::ostringstream ss;
    ss << rootPath << "/" << iterationIdx << "/";

    return ss.str();
}


void OptimizationResultWriter::logParPEVersion() {
    hdf5WriteStringAttribute(file_id, rootPath.c_str(), "PARPE_VERSION",
                             GIT_VERSION);
}


void OptimizationResultWriter::closeResultHDFFile() {
    H5_SAVE_ERROR_HANDLER;
    herr_t status = H5Fclose(file_id);

    if (status < 0) {
        error("closeResultHDFFile failed to close HDF5 file.");
        H5Eprint(H5E_DEFAULT, stderr);
    }
    H5_RESTORE_ERROR_HANDLER;
}


void OptimizationResultWriter::logLocalOptimizerObjectiveFunctionEvaluation(
        const double *parameters, int numParameters, double objectiveFunctionValue,
        const double *objectiveFunctionGradient, int numIterations, int numFunctionCalls,
        double timeElapsedInSeconds) {

    std::string pathStr = getIterationPath(numIterations);
    const char *fullGroupPath = pathStr.c_str();

    hdf5CreateOrExtendAndWriteToDouble2DArray(
                file_id, fullGroupPath, "costFunCost", &objectiveFunctionValue, 1);

    if (objectiveFunctionGradient && numParameters) {
        hdf5CreateOrExtendAndWriteToDouble2DArray(
                    file_id, fullGroupPath, "costFunGradient",
                    objectiveFunctionGradient, numParameters);
    } else if (numParameters) {
        double dummyGradient[numParameters];
        std::fill_n(dummyGradient, numParameters, NAN);
        hdf5CreateOrExtendAndWriteToDouble2DArray(file_id, fullGroupPath,
                                                  "costFunGradient",
                                                  dummyGradient, numParameters);
    }

    if (parameters && numParameters)
        hdf5CreateOrExtendAndWriteToDouble2DArray(file_id, fullGroupPath,
                                                  "costFunParameters",
                                                  parameters, numParameters);

    hdf5CreateOrExtendAndWriteToDouble2DArray(file_id, fullGroupPath,
                                              "costFunWallTimeInSec",
                                              &timeElapsedInSeconds, 1);

    hdf5CreateOrExtendAndWriteToInt2DArray(
                file_id, fullGroupPath, "costFunCallIndex", &numFunctionCalls, 1);

    flushResultWriter();
}


void OptimizationResultWriter::logLocalOptimizerIteration(
        int numIterations, const double * const theta, int numParameters,
        double objectiveFunctionValue, const double * const gradient,
        double timeElapsedInSeconds) {
    std::string pathStr = getOptimizationPath();
    const char *fullGroupPath = pathStr.c_str();

    hdf5CreateOrExtendAndWriteToDouble2DArray(
                file_id, fullGroupPath, "iterCostFunCost", &objectiveFunctionValue, 1);
    if (gradient)
        hdf5CreateOrExtendAndWriteToDouble2DArray(file_id, fullGroupPath,
                                                  "iterCostFunGradient",
                                                  gradient, numParameters);
    if (theta)
        hdf5CreateOrExtendAndWriteToDouble2DArray(file_id, fullGroupPath,
                                                  "iterCostFunParameters",
                                                  theta, numParameters);
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


void OptimizationResultWriter::starting(int numParameters, const double * const initialParameters)
{
    if (initialParameters) {
        std::string pathStr = getOptimizationPath();
        const char *fullGroupPath = pathStr.c_str();

        hdf5CreateOrExtendAndWriteToDouble2DArray(file_id, fullGroupPath,
                                                  "initialParameters",
                                                  initialParameters, numParameters);
    }
}


void OptimizationResultWriter::flushResultWriter() {
    auto lock = hdf5MutexGetLock();

    H5Fflush(file_id, H5F_SCOPE_LOCAL);
}

void OptimizationResultWriter::saveLocalOptimizerResults(
        double finalNegLogLikelihood, const double *optimalParameters,
        int numParameters, double masterTime, int exitStatus) {

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

    if (optimalParameters) {
        fullGroupPath = (optimPath + "/finalParameters");
        dimensions[0] = numParameters;
        H5LTmake_dataset(file_id, fullGroupPath.c_str(), 1, dimensions,
                         H5T_NATIVE_DOUBLE, optimalParameters);
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
