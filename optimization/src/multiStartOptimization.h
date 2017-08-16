#ifndef MULTI_START_OPTIMIZATION_H
#define MULTI_START_OPTIMIZATION_H

#include "optimizationProblem.h"

class OptimizationProblemGeneratorForMultiStart {

  public:
    OptimizationProblemGeneratorForMultiStart();

    OptimizationProblem *getLocalProblem(int multiStartIndex);

    OptimizationProblem **
    createLocalOptimizationProblems(int numLocalOptimizations);

    virtual OptimizationProblem *getLocalProblemImpl(int multiStartIndex) = 0;
};

/**
 * @brief runParallelMultiStartOptimization
 * @return always returns 0
 */

int runParallelMultiStartOptimization(
    OptimizationProblemGeneratorForMultiStart *problemGenerator,
    int numberOfStarts, bool restartOnFailure);

#endif
