#include "optimizationResultWriter.h"

#include <parpeVersion.h>
#include <hdf5Misc.h>
#include <misc.h>
#include <logging.h>

#include <cmath>
#include <sstream>
#include <iostream>

#include <hdf5_hl.h>


namespace parpe {

OptimizationResultWriter::OptimizationResultWriter(hid_t file_id, std::string rootPath)
    : rootPath(rootPath)
{
    auto lock = hdf5MutexGetLock();
    this->file_id = H5Freopen(file_id);

    hdf5EnsureGroupExists(file_id, this->rootPath);
}


OptimizationResultWriter::OptimizationResultWriter(const std::string &filename,
                                                   bool overwrite, std::string rootPath)
    : rootPath(std::move(rootPath))
{
    logmessage(LOGLVL_DEBUG, "Writing results to %s.", filename.c_str());

    file_id = hdf5CreateFile(filename.c_str(), overwrite);

    hdf5EnsureGroupExists(file_id, this->rootPath);
}

OptimizationResultWriter::OptimizationResultWriter(const OptimizationResultWriter &other)
    : rootPath(other.rootPath)
{
    auto lock = hdf5MutexGetLock();
    file_id = H5Freopen(other.file_id);
    hdf5EnsureGroupExists(file_id, rootPath);
}


const std::string &OptimizationResultWriter::getRootPath() const { return rootPath; }


std::string OptimizationResultWriter::getIterationPath(int iterationIdx) const {
    std::ostringstream ss;
    ss << rootPath << "/iteration/" << iterationIdx << "/";

    return ss.str();
}


void OptimizationResultWriter::logObjectiveFunctionEvaluation(
        gsl::span<const double> parameters,
        double objectiveFunctionValue,
        gsl::span<const double> objectiveFunctionGradient,
        int numIterations, int numFunctionCalls,
        double timeElapsedInSeconds, 
		bool logGradientEachFunctionEvaluation,
		bool logParametersEachFunctionEvaluation)
{

    std::string pathStr = getIterationPath(numIterations);
    const char *fullGroupPath = pathStr.c_str();

    hdf5CreateOrExtendAndWriteToDouble2DArray(
                file_id, fullGroupPath, "costFunCost", &objectiveFunctionValue, 1);

    if (logGradientEachFunctionEvaluation) {
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
    }

    if (logParametersEachFunctionEvaluation)
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


void OptimizationResultWriter::logOptimizerIteration(int numIterations,
                                                     gsl::span<const double> parameters,
                                                     double objectiveFunctionValue,
                                                     gsl::span<const double> gradient,
                                                     double wallSeconds,
                                                     double cpuSeconds,
													 bool logGradientEachIteration)
{
    std::string const& pathStr = getRootPath();
    const char *fullGroupPath = pathStr.c_str();

    hdf5CreateOrExtendAndWriteToDouble2DArray(
                file_id, fullGroupPath, "iterCostFunCost", &objectiveFunctionValue, 1);

    if (logGradientEachIteration) {
		if (!gradient.empty()) {
			hdf5CreateOrExtendAndWriteToDouble2DArray(file_id, fullGroupPath,
													  "iterCostFunGradient",
													  gradient.data(), gradient.size());
		} else if(!parameters.empty()) {
			std::vector<double> nanGradient(parameters.size(), NAN);
			hdf5CreateOrExtendAndWriteToDouble2DArray(file_id, fullGroupPath,
													  "iterCostFunGradient",
													  nanGradient.data(), nanGradient.size());
		}
    }

    if (!parameters.empty()) {
        hdf5CreateOrExtendAndWriteToDouble2DArray(file_id, fullGroupPath,
                                                  "iterCostFunParameters",
                                                  parameters.data(), parameters.size());
    }

    hdf5CreateOrExtendAndWriteToDouble2DArray(file_id, fullGroupPath,
                                              "iterCostFunWallSec",
                                              &wallSeconds, 1);

    hdf5CreateOrExtendAndWriteToDouble2DArray(file_id, fullGroupPath,
                                              "iterCostFunCpuSec",
                                              &cpuSeconds, 1);

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
        std::string const& pathStr = getRootPath();
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

void OptimizationResultWriter::saveOptimizerResults(double finalNegLogLikelihood,
                                                    gsl::span<const double> optimalParameters,
                                                    double wallSec,
                                                    double cpuSec,
                                                    int exitStatus) const
{

    std::string const& optimPath = getRootPath();
    hdf5EnsureGroupExists(file_id, optimPath.c_str());

    std::string fullGroupPath;
    hsize_t dimensions[1] = {1};

    auto lock = hdf5MutexGetLock();

    fullGroupPath = (optimPath + "/finalCost");
    H5LTmake_dataset(file_id, fullGroupPath.c_str(), 1, dimensions,
                     H5T_NATIVE_DOUBLE, &finalNegLogLikelihood);

    fullGroupPath = (optimPath + "/wallSec");
    H5LTmake_dataset(file_id, fullGroupPath.c_str(), 1, dimensions,
                     H5T_NATIVE_DOUBLE, &wallSec);

    fullGroupPath = (optimPath + "/cpuSec");
    H5LTmake_dataset(file_id, fullGroupPath.c_str(), 1, dimensions,
                     H5T_NATIVE_DOUBLE, &cpuSec);

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
    closeHDF5File(file_id);
}

hid_t OptimizationResultWriter::getFileId() const
{
    return file_id;
}

} // namespace parpe
