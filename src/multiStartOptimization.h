#ifndef MULTI_START_OPTIMIZATION_H
#define MULTI_START_OPTIMIZATION_H

#include "optimizationProblem.h"

typedef struct MultiStartOptimization_tag MultiStartOptimization;

typedef int (*multiStartOptimizationGetInitialPointFp)(const MultiStartOptimization *multiStartOptimization, int currentStartIdx, double *buffer);
typedef void *(*multiStartOptimizationGetUserDataFp)(const MultiStartOptimization *multiStartOptimization, int currentStartIdx);


// TODO: add userData for multistart IDs, ...
typedef struct MultiStartOptimization_tag {

    OptimizationProblem *optimizationProblem;

    /** Number of starts to perform */
    int numberOfStarts;

    /** restart local optimization, until numberOfStarts successful starts are reached */
    bool restartOnFailure;

    /** Function for initial points. If 0, random starting points are drawn from [optimizationProblem.parametersMin, optimizationProblem.parametersMax] */
    multiStartOptimizationGetInitialPointFp getInitialPoint;

    /** Function to get userData for each start */
    multiStartOptimizationGetUserDataFp getUserData;
} MultiStartOptimization;


MultiStartOptimization *multiStartOptimizationNew();

/**
 * @brief runParallelMultiStartOptimization
 * @param multiStartOptimization
 * @return always returns 0
 */

int runParallelMultiStartOptimization(MultiStartOptimization *multiStartOptimization);

#endif
