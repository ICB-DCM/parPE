#ifndef MULTI_START_OPTIMIZATION_H
#define MULTI_START_OPTIMIZATION_H

#include <parpeoptimization/optimizationProblem.h>

#include <memory>

namespace parpe {

/**
 * @brief Interface for multi-start optimization problems
 */
class MultiStartOptimizationProblem {
public:
    virtual int getNumberOfStarts() const = 0;

    virtual bool restartOnFailure() const { return false; }

    virtual std::unique_ptr<OptimizationProblem>
    getLocalProblem(int multiStartIndex) const = 0;

    virtual ~MultiStartOptimizationProblem() = default;
};


/**
 * @brief The MultiStartOptimization class runs multiple optimization runs
 */

class MultiStartOptimization {

  public:
    MultiStartOptimization(MultiStartOptimizationProblem& problem,
                           bool runParallel = true,
                           int first_start_idx = 0);

    ~MultiStartOptimization() = default;

    /**
     * @brief Start multi-start optimization
     */
    void run();

    /**
     * @brief Run all optimizations in parallel, each in a dedicated thread
     */
    void runMultiThreaded() const;

    /**
     * @brief Run optimizations sequentially
     */
    void runSingleThreaded();

    /**
     * @brief Set parallel or sequential mode
     */
    void setRunParallel(bool runParallel);

  private:
    /**
     * @brief Optimize local problem for the given start index
     */
    int runStart(int start_idx) const;


    /** Optimization problem to be solved */
    MultiStartOptimizationProblem& msProblem;

    /** Number of optimization runs to perform*/
    int numberOfStarts = 1;

    /** Try a new starting point if a previous one fails */
    bool restartOnFailure = false;

    /** Run multiple optimizations in parallel */
    bool runParallel = true;

    /** Index value of the first start
     * Usable when splitting starts across multiple files */
    int first_start_idx = 0;
};

} // namespace parpe

#endif
