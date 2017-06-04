#include "multiConditionProblemResultWriter.h"
#include <assert.h>
#include <logging.h>
#include "hdf5Misc.h"
#include <sys/stat.h>

#define HDF5_PATH_BUFFER_MAX_LEN 1024

// HDF5 file to write results to
static hid_t file_id = 0;

int initResultHDFFile(const char *filename, bool overwrite)
{
    H5_SAVE_ERROR_HANDLER;;

    if(!overwrite) {
        struct stat st = {0};
        bool fileExists = stat(filename, &st) == -1;

        assert(fileExists);
    }

    // overwrites TODO: backup
    file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    if(file_id < 0) {
        H5Eprint(H5E_DEFAULT, stderr);
    }

    H5_RESTORE_ERROR_HANDLER;;

    return file_id < 0;
}


void closeResultHDFFile()
{
    H5_SAVE_ERROR_HANDLER;
    herr_t status = H5Fclose(file_id);

    if(status< 0) {
        error("closeResultHDFFile failed to close HDF5 file.");
        H5Eprint(H5E_DEFAULT, stderr);
    }
    H5_RESTORE_ERROR_HANDLER;
}


void flushResultWriter()
{
    hdf5LockMutex();

    H5Fflush(file_id, H5F_SCOPE_LOCAL);

    hdf5UnlockMutex();
}



void logLocalOptimizerObjectiveFunctionEvaluation(OptimizationProblem *problem,
                                                  const double *parameters,
                                                  double objectiveFunctionValue,
                                                  const double *gradient,
                                                  int numFunctionCalls,
                                                  double timeElapsedInSeconds)
{
    char fullGroupPath[HDF5_PATH_BUFFER_MAX_LEN];
    sprintf(fullGroupPath, "/multistarts/%d/objFunGradEval/",
            1); // TODO multistartindex idxLocalOptimization

    hdf5CreateOrExtendAndWriteToDouble2DArray(file_id, fullGroupPath, "negLogLikelihood", &objectiveFunctionValue, 1);
    if(gradient)
        hdf5CreateOrExtendAndWriteToDouble2DArray(file_id, fullGroupPath, "negLogLikelihoodGradient", gradient, problem->numOptimizationParameters);
    hdf5CreateOrExtendAndWriteToDouble2DArray(file_id, fullGroupPath, "p", parameters, problem->numOptimizationParameters);
    hdf5CreateOrExtendAndWriteToDouble2DArray(file_id, fullGroupPath, "evalFTime", &timeElapsedInSeconds, 1);

    flushResultWriter();
}

void logLocalOptimizerIteration(JobIdentifier id, int numIterations, double *theta, double objectiveFunctionValue, const double *gradient, double timeElapsedInSeconds, int nTheta, int alg_mod, double inf_pr, double inf_du, double mu, double d_norm, double regularization_size, double alpha_du, double alpha_pr, int ls_trials)
{
    char fullGroupPath[HDF5_PATH_BUFFER_MAX_LEN];
    sprintf(fullGroupPath, "/multistarts/%d/optimizerIterations/", id.idxLocalOptimization);

    hdf5CreateOrExtendAndWriteToDouble2DArray(file_id, fullGroupPath, "negLogLikelihood", &objectiveFunctionValue, 1);
    hdf5CreateOrExtendAndWriteToDouble2DArray(file_id, fullGroupPath, "negLogLikelihoodGradient", gradient, nTheta);
    hdf5CreateOrExtendAndWriteToDouble2DArray(file_id, fullGroupPath, "p", theta, nTheta);
    hdf5CreateOrExtendAndWriteToDouble2DArray(file_id, fullGroupPath, "evalFTime", &timeElapsedInSeconds, 1);

    hdf5CreateOrExtendAndWriteToInt2DArray(   file_id, fullGroupPath, "numIterations", &numIterations, 1);
    hdf5CreateOrExtendAndWriteToInt2DArray(   file_id, fullGroupPath, "alg_mod", &alg_mod, 1);
    hdf5CreateOrExtendAndWriteToDouble2DArray(file_id, fullGroupPath, "inf_pr", &inf_pr, 1);
    hdf5CreateOrExtendAndWriteToDouble2DArray(file_id, fullGroupPath, "inf_du", &inf_du, 1);
    hdf5CreateOrExtendAndWriteToDouble2DArray(file_id, fullGroupPath, "mu", &mu, 1);
    hdf5CreateOrExtendAndWriteToDouble2DArray(file_id, fullGroupPath, "d_norm", &d_norm, 1);
    hdf5CreateOrExtendAndWriteToDouble2DArray(file_id, fullGroupPath, "regularization_size", &regularization_size, 1);
    hdf5CreateOrExtendAndWriteToDouble2DArray(file_id, fullGroupPath, "alpha_du", &alpha_du, 1);
    hdf5CreateOrExtendAndWriteToDouble2DArray(file_id, fullGroupPath, "alpha_pr", &alpha_pr, 1);
    hdf5CreateOrExtendAndWriteToInt2DArray(   file_id, fullGroupPath, "ls_trials", &ls_trials, 1);

    flushResultWriter();
}

void logSimulation(JobIdentifier id, const double *theta, double llh, const double *gradient, double timeElapsedInSeconds, int nTheta, int numStates, double *states,
                   double *stateSensi, double *y, int jobId, int iterationsUntilSteadystate)
{
    char fullGroupPath[HDF5_PATH_BUFFER_MAX_LEN];
    // TODO doesnt have to be extendable in current implementation
    sprintf(fullGroupPath, "/multistarts/%d/iteration/%d/genotype/%d/experiment/%d/",
            id.idxLocalOptimization, id.idxLocalOptimizationIteration, 0, id.idxConditions);

    hdf5CreateOrExtendAndWriteToDouble2DArray(file_id, fullGroupPath, "negLogLikelihood", &llh, 1);
    hdf5CreateOrExtendAndWriteToInt2DArray(file_id, fullGroupPath, "jobId", &jobId, 1);

    if(gradient)
        hdf5CreateOrExtendAndWriteToDouble2DArray(file_id, fullGroupPath, "negLogLikelihoodGradient", gradient, nTheta);

    if(theta)
        hdf5CreateOrExtendAndWriteToDouble2DArray(file_id, fullGroupPath, "p", theta, nTheta);
    hdf5CreateOrExtendAndWriteToDouble2DArray(file_id, fullGroupPath, "evalFTime", &timeElapsedInSeconds, 1);

    if(states)
        hdf5CreateOrExtendAndWriteToDouble2DArray(file_id, fullGroupPath, "X", states, numStates);
    if(y)
        hdf5CreateOrExtendAndWriteToDouble2DArray(file_id, fullGroupPath, "Y", y, 1);

    hdf5CreateOrExtendAndWriteToInt2DArray(   file_id, fullGroupPath, "iterationsUntilSteadystate", &iterationsUntilSteadystate, 1);

    if(stateSensi)
        hdf5CreateOrExtendAndWriteToDouble3DArray(file_id, fullGroupPath, "sX", stateSensi, numStates, nTheta);

    flushResultWriter();
}


void saveTotalWalltime(const double timeInSeconds) {
    hsize_t dims[1] = {1};

    hdf5LockMutex();

    H5LTmake_dataset(file_id, "/totalWallTimeInSec", 1, dims, H5T_NATIVE_DOUBLE, &timeInSeconds);

    hdf5UnlockMutex();
}

void saveLocalOptimizerResults(OptimizationProblem *problem,  JobIdentifier path,
                               double finalNegLogLikelihood,
                               const double *optimalParameters,
                               double timeElapsed,
                               int exitStatus)
{
    hsize_t dimensions[1] = {1};

    char fullPath[HDF5_PATH_BUFFER_MAX_LEN];

//    hdf5LockMutex();

//    sprintf(fullPath, "/crossvalidations/%d/multistarts/%d/finalNegLogLikelihood", path.idxMultiStart, path.idxLocalOptimization);
//    H5LTmake_dataset(file_id, fullPath, 1, dimensions, H5T_NATIVE_DOUBLE, &finalNegLogLikelihood);

//    sprintf(fullPath, "/crossvalidations/%d/multistarts/%d/timeOnMaster", path.idxMultiStart, path.idxLocalOptimization);
//    H5LTmake_dataset(file_id, fullPath, 1, dimensions, H5T_NATIVE_DOUBLE, &timeElapsed);

//    sprintf(fullPath, "/crossvalidations/%d/multistarts/%d/exitStatus", path.idxMultiStart, path.idxLocalOptimization);
//    H5LTmake_dataset(file_id, fullPath, 1, dimensions, H5T_NATIVE_INT, &exitStatus);

//    dimensions[0] = problem->numOptimizationParameters;
//    sprintf(fullPath, "/crossvalidations/%d/multistarts/%d/optimalParameters", path.idxMultiStart, path.idxLocalOptimization);
//    H5LTmake_dataset(file_id, fullPath, 1, dimensions, H5T_NATIVE_DOUBLE, &optimalParameters);

//    hdf5UnlockMutex();

    flushResultWriter();
}
