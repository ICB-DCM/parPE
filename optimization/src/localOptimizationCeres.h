#ifndef LOCAL_OPTIMIZATION_CERES_H
#define LOCAL_OPTIMIZATION_CERES_H

#include "optimizationProblem.h"

#ifdef __cplusplus
#define EXTERNC extern "C"
#else
#define EXTERNC
#endif

/**
 * @brief getLocalOptimumCeres determines the local optimum for the provided
 * optimization problem using the Google Ceres optimizer
 * @param problem the optimization problem
 * @return Returns 0 on success.
 */

EXTERNC int getLocalOptimumCeres(OptimizationProblem *problem);

#endif
