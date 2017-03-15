#ifndef LOGGER_H
#define LOGGER_H

#include <hdf5.h>
#include <hdf5_hl.h>
#include <stdbool.h>
#include "../objectiveFunctionBenchmarkModel/dataprovider.h"

int initResultHDFFile(const char *filename);

void closeResultHDFFile();

void logLocalOptimizerObjectiveFunctionEvaluation(datapath path, int numFunctionCalls, double *theta, double objectiveFunctionValue, double timeElapsedInSeconds, int nTheta);
void logLocalOptimizerObjectiveFunctionGradientEvaluation(datapath path, int numFunctionCalls, double *theta, double objectiveFunctionValue, const double *gradient, double timeElapsedInSeconds, int nTheta);

// TODO add likelihood for individual experiments
void logLocalOptimizerIteration(datapath path, int numIterations, double *theta, double objectiveFunctionValue, const double *gradient, double timeElapsedInSeconds, int nTheta, int alg_mod, double inf_pr, double inf_du, double mu, double d_norm, double regularization_size, double alpha_du, double alpha_pr, int ls_trials);

void logSimulation(datapath path, double llh, const double *gradient, double timeElapsedInSeconds, int nTheta, int numStates, double *states, double *stateSensi, double *y, int jobId, int iterationsUntilSteadystate);

void flushResultWriter();

void saveTotalWalltime(const double timeInSeconds);

void saveLocalOptimizerResults(datapath path, double finalNegLogLikelihood, double masterTime, int exitStatus);

#endif
