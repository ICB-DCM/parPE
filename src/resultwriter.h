#ifndef LOGGER_H
#define LOGGER_H

#include <hdf5.h>
#include <hdf5_hl.h>
#include <stdbool.h>
#include "dataprovider.h"

#define LOGGER_H_LOGFILE "mylog.log"

typedef struct {
    hid_t file_id;
    char *rootPath;

} loggerdata;


int initResultHDFFile(const char *filename);

void closeResultHDFFile();

void logLocalOptimizerObjectiveFunctionEvaluation(datapath path, int numFunctionCalls, double *theta, double objectiveFunctionValue, double timeElapsedInSeconds, int nTheta);
void logLocalOptimizerObjectiveFunctionGradientEvaluation(datapath path, int numFunctionCalls, double *theta, double objectiveFunctionValue, const double *gradient, double timeElapsedInSeconds, int nTheta);

// TODO add likelihood for individual experiments
void logLocalOptimizerIteration(datapath path, int numIterations, double *theta, double objectiveFunctionValue, const double *gradient, double timeElapsedInSeconds, int nTheta, int alg_mod, double inf_pr, double inf_du, double mu, double d_norm, double regularization_size, double alpha_du, double alpha_pr, int ls_trials);

void logSimulation(datapath path, double llh, const double *gradient, double timeElapsedInSeconds, int nTheta, int numStates, double *states, double *stateSensi, double *y, int jobId, int iterationsUntilSteadystate);

void flushResultWriter();

// TODO
bool hdf5DatasetExists(hid_t file_id, const char *datasetName);

bool hdf5GroupExists(hid_t file_id, const char *groupName);

void hdf5CreateGroup(hid_t file_id, const char *groupPath, bool recursively);

void hdf5CreateExtendableDouble2DArray(hid_t file_id, const char* datasetPath, int stride);

void hdf5CreateExtendableInt2DArray(hid_t file_id, const char* datasetPath, int stride);

void hdf5CreateExtendableDouble3DArray(hid_t file_id, const char *datasetPath, int stride1, int stride2);

void hdf5Extend2ndDimensionAndWriteToDouble2DArray(hid_t file_id, const char *datasetPath, const double *buffer);

void hdf5Extend2ndDimensionAndWriteToInt2DArray(hid_t file_id, const char *datasetPath, const int *buffer);

void hdf5CreateOrExtendAndWriteToDouble2DArray(hid_t file_id, const char *parentPath, const char *datasetName, const double *buffer, int stride);

void hdf5CreateOrExtendAndWriteToInt2DArray(hid_t file_id, const char *parentPath, const char *datasetName, const int *buffer, int stride);

void hdf5CreateOrExtendAndWriteToDouble3DArray(hid_t file_id, const char *parentPath, const char *datasetName, const double *buffer, int stride1, int stride2);

char *myStringCat(const char *first, const char *second);

void saveTotalWalltime(const double timeInSeconds);

void saveLocalOptimizerResults(datapath path, double finalNegLogLikelihood, double masterTime, int exitStatus);

#endif
