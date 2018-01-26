#ifndef MULTI_START_OPTIMIZATION_H
#define MULTI_START_OPTIMIZATION_H

#include "optimizationProblem.h"
#include <vector>
#include <memory>

namespace parpe {

/**
 * @brief Interface for multi-start optimization problems
 */
class MultiStartOptimizationProblem {
public:
    virtual int getNumberOfStarts() const = 0;

    virtual bool restartOnFailure() const { return false; }

    virtual std::unique_ptr<OptimizationProblem> getLocalProblem(int multiStartIndex) const = 0;
};


/**
 * @brief The MultiStartOptimization class runs multiple optimization runs
 */

class MultiStartOptimization {

  public:
    MultiStartOptimization(MultiStartOptimizationProblem& problem, bool runParallel = true);

    /**
     * @brief Start multi-start optimization
     * @return always returns 0
     */

    int run();

    int runMultiThreaded();

    int runSingleThreaded();

    void setRunParallel(bool runParallel);

  private:

    std::vector<OptimizationProblem *> createLocalOptimizationProblems();

    MultiStartOptimizationProblem& msProblem;
    int numberOfStarts = 1;
    bool restartOnFailure = false;
    bool runParallel = true;
};

} // namespace parpe

#endif
