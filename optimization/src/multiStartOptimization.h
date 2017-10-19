#ifndef MULTI_START_OPTIMIZATION_H
#define MULTI_START_OPTIMIZATION_H

#include "optimizationProblem.h"
#include <vector>

namespace parpe {

class MultiStartOptimization {

  public:
    MultiStartOptimization() = default;

    MultiStartOptimization(int numberOfStarts, bool restartOnFailure);

    /**
     * @brief Start multi-start optimization
     * @return always returns 0
     */

    int run();

  protected:
    OptimizationProblem *getLocalProblem(int multiStartIndex);

    std::vector<OptimizationProblem *> createLocalOptimizationProblems();

    virtual OptimizationProblem *getLocalProblemImpl(int multiStartIndex) = 0;

    int numberOfStarts = 1;
    bool restartOnFailure = false;
};

} // namespace parpe

#endif
