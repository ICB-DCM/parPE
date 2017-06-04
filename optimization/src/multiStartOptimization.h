#ifndef MULTI_START_OPTIMIZATION_H
#define MULTI_START_OPTIMIZATION_H

#include "optimizationProblem.h"

typedef OptimizationProblem* (*optimizationProblemGeneratorForMultiStartFp)(int currentStartIdx, void *userData);

/**
 * @brief runParallelMultiStartOptimization
 * @return always returns 0
 */

int runParallelMultiStartOptimization(optimizationProblemGeneratorForMultiStartFp problemGenerator,
                                      int numberOfStarts, bool restartOnFailure, void *userData);

#endif
