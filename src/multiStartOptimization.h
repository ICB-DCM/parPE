#ifndef MULTI_START_OPTIMIZATION_H
#define MULTI_START_OPTIMIZATION_H

#include "optimizationProblem.h"

typedef int (*multiStartOptimizationGetInitialPointFp)(void *multiStartOptimization, int currentStartIdx, double *buffer);
typedef void *(*multiStartOptimizationGetUserDataFp)(void *multiStartOptimization, int currentStartIdx);


// TODO: add userData for multistart IDs, ...
typedef struct MultiStartOptimization_tag {

    OptimizationProblem *optimizationProblem;

    /** Number of starts to perform */
    int numberOfStarts;

    /** restart local optimization, until numberOfStarts successful starts are reached */
    bool restartOnFailure;

    /** Function for initial points */
    multiStartOptimizationGetInitialPointFp getInitialPoint;

    /** Function to get userData for each start */
    multiStartOptimizationGetUserDataFp getUserData;
} MultiStartOptimization;


MultiStartOptimization *multiStartOptimizationNew();

int runParallelMultiStartOptimization(MultiStartOptimization *multiStartOptimization);


#endif
