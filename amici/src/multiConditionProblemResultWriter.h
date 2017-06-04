#ifndef LOGGER_H
#define LOGGER_H

#include <hdf5_hl.h>
#include <stdbool.h>
#include "optimizationProblem.h"
#include "MultiConditionDataProvider.h"

int initResultHDFFile(const char *filename, bool overwrite);

void closeResultHDFFile();

void logLocalOptimizerObjectiveFunctionEvaluation(OptimizationProblem *problem,
                                                  const double *parameters,
                                                  double objectiveFunctionValue,
                                                  const double *objectiveFunctionGradient,
                                                  int numFunctionCalls,
                                                  double timeElapsed);

// TODO add likelihood for individual experiments
void logLocalOptimizerIteration(JobIdentifier path,
                                int numIterations,
                                double *theta,
                                double objectiveFunctionValue,
                                const double *gradient,
                                double timeElapsedInSeconds,
                                int nTheta,
                                int alg_mod,
                                double inf_pr, double inf_du,
                                double mu,
                                double d_norm,
                                double regularization_size,
                                double alpha_du, double alpha_pr,
                                int ls_trials);

void logSimulation(JobIdentifier path, const double *theta,
                   double llh, const double *gradient,
                   double timeElapsedInSeconds,
                   int nTheta, int numStates,
                   double *states, double *stateSensi,
                   double *y, int jobId,
                   int iterationsUntilSteadystate);

void flushResultWriter();

void saveTotalWalltime(const double timeInSeconds);

void saveLocalOptimizerResults(OptimizationProblem *problem, JobIdentifier path, double finalNegLogLikelihood,
                               const double *optimalParameters,
                               double masterTime,
                               int exitStatus);

#endif
