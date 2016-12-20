#ifndef DATA_PROVIDER_H
#define DATA_PROVIDER_H

#define EXPERIMENT_INDEX_CONTROL -1
#define XDOT_REL_TOLERANCE 1e-6
#define NUM_OPTIMIZATION_PARAMS 4141
#define NUM_FIXED_PARAMS 765
#define NUM_STATE_VARIABLES 1230

#define DATA_PROVIDER_H_VERBOSE 0

#include <hdf5.h>
#include <hdf5_hl.h>
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include "misc.h"

typedef struct datapath_tag {
    int idxMultiStart; // current multistart batch (e.g. for crossvalidation)
    int idxLocalOptimization; // current start index in multistart run
    int idxLocalOptimizationIteration; // iteration of local solver
    int idxGenotype;
    int idxExperiment;
} datapath;

void printDatapath(datapath p);

void sprintDatapath(char *buffer, datapath path);

int initDataProvider(const char *hdf5Filename);

void closeDataProvider();

int getNumMultiStartRuns();

int getNumLocalOptimizationsForMultiStartRun(int multiStartRun);

int getMaxIter();

int getNumGenotypes(datapath path);

int getExperimentCountForCellline(datapath dpath);

int getLenTheta();

int getLenKappa();

void getThetaLowerBounds(datapath dataPath, double *buffer, AMI_parameter_scaling scaling);

void getThetaUpperBounds(datapath dataPath, double *buffer, AMI_parameter_scaling scaling);

void getInitialTheta(datapath dataPath, double *buffer, AMI_parameter_scaling scaling);

void getFeasibleInitialThetaFromFile(datapath dataPath, double *buffer, AMI_parameter_scaling scaling);

void getRandomInitialThetaFromFile(datapath dataPath, double *buffer, AMI_parameter_scaling scaling);

int getNumberOfSimulationForObjectiveFunction(datapath path);

UserData *getMyUserData(); // TODO: specific for each multistart run

void readFixedParameters(int dataMatrixRow, UserData *udata);

ExpData *getExperimentalDataForExperiment(datapath dpath, UserData *udata);

// TODO void logmessage(loglevel lvl, datapath path, const char *format, ...);


// Workaround to delete ExpData memory allocated by getExperimentalDataForExperiment, since AMICI freeExpData cannot be used (new[] vs malloc issue)
void myFreeExpData(ExpData *edata);

void freeUserDataC(UserData *udata);

herr_t hdf5ErrorStackWalker_cb(unsigned int n, const H5E_error_t *err_desc, void *client_data);

// malloc version of ami_hdf5.cpp
void getDoubleArrayAttributeC(hid_t file_id, const char* optionsObject, const char* attributeName, double **destination, hsize_t *length);

void getIntArrayAttributeC(hid_t file_id, const char* optionsObject, const char* attributeName, int **destination, hsize_t *length);

void dataproviderPrintInfo();

#endif
