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

    virtual ~MultiStartOptimizationProblem() = default;
};


/**
 * @brief The MultiStartOptimization class runs multiple optimization runs
 */

class MultiStartOptimization {

  public:
    MultiStartOptimization(MultiStartOptimizationProblem& problem,
                           bool runParallel = true);

    ~MultiStartOptimization() = default;

    /**
     * @brief Start multi-start optimization
     */
    void run();

    /**
     * @brief Run all optimizations in parallel, each in a dedicated thread
     */
    void runMultiThreaded();

    /**
     * @brief Run optimizations sequentially
     */
    void runSingleThreaded();

    /**
     * @brief Set parallel or sequential mode
     */
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
