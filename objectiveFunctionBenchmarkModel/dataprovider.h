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
#include <include/udata.h>
#include <include/edata.h>
#include "misc.h"

typedef struct datapath_tag {
    int idxMultiStart;                 /** current multistart batch (e.g. for crossvalidation) */
    int idxLocalOptimization;          /** current start index in multistart run */
    int idxLocalOptimizationIteration; /** iteration of local solver */
    int idxGenotype;                   /** genotype index (starting from 1!!) */
    int idxExperiment;                 /** experiment index */
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

/**
 * @brief getThetaLowerBounds Get lower parameter bounds
 * @param dataPath (not yet used)
 * @param buffer allocated memory to write parameter bounds
 * @param scaling scaling methods for parameter bounds
 */

void getThetaLowerBounds(datapath dataPath, double *buffer, AMI_parameter_scaling scaling);

/**
 * @brief getThetaUpperBounds Get upper parameter bounds
 * @param dataPath (not yet used)
 * @param buffer allocated memory to write parameter bounds
 * @param scaling scaling methods for parameter bounds
 */

void getThetaUpperBounds(datapath dataPath, double *buffer, AMI_parameter_scaling scaling);

/**
 * @brief getInitialThetaLHS Get random starting parameters using latin hypercube sampling. TODO: check LHS code
 * @param dataPath
 * @param buffer
 * @param scaling
 */

void getInitialThetaLHS(datapath dataPath, double *buffer, AMI_parameter_scaling scaling);

/**
 * @brief getFeasibleInitialThetaFromFile Get user-provided initial theta TODO: can be removed?
 * @param dataPath
 * @param buffer
 * @param scaling
 */

void getFeasibleInitialThetaFromFile(datapath dataPath, double *buffer, AMI_parameter_scaling scaling);

/**
 * @brief getRandomInitialThetaFromFile Get a initial parameter vector from a pregenerated list of random vectors
 * @param dataPath
 * @param buffer
 * @param scaling NOT USED
 */

void getRandomInitialThetaFromFile(datapath dataPath, double *buffer, AMI_parameter_scaling scaling);

int getNumberOfSimulationForObjectiveFunction(datapath path); // Not used

UserData *getMyUserData(); // TODO: specific for each multistart run

ExpData *getExperimentalDataForExperiment(datapath dpath, UserData *udata);

herr_t hdf5ErrorStackWalker_cb(unsigned int n, const H5E_error_t *err_desc, void *client_data); // TODO: also use for resultwriter

void dataproviderPrintInfo();

// TODO void logmessage(loglevel lvl, datapath path, const char *format, ...);

// Workaround to delete ExpData memory allocated by getExperimentalDataForExperiment, since AMICI freeExpData cannot be used (new[] vs malloc issue)
void myFreeExpData(ExpData *edata);

void freeUserDataC(UserData *udata);

#endif
