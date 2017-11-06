#ifndef MULTI_START_OPTIMIZATION_H
#define MULTI_START_OPTIMIZATION_H

#include "optimizationProblem.h"
#include <vector>
#include <memory>

namespace parpe {

class MultiStartOptimization {

  public:
    MultiStartOptimization() = default;

    MultiStartOptimization(int numberOfStarts, bool restartOnFailure, bool runParallel = true);

    MultiStartOptimization(OptimizationOptions const& o);

    /**
     * @brief Start multi-start optimization
     * @return always returns 0
     */

    int run();

    int runMultiThreaded();

    int runSingleThreaded();

    void setRunParallel(bool runParallel);

  protected:

    std::unique_ptr<OptimizationProblem> getLocalProblem(int multiStartIndex);

    std::vector<OptimizationProblem *> createLocalOptimizationProblems();

    virtual std::unique_ptr<OptimizationProblem> getLocalProblemImpl(int multiStartIndex) = 0;

    int numberOfStarts = 1;
    bool restartOnFailure = false;
    bool runParallel = true;
};

} // namespace parpe

#endif
