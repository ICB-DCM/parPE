#ifndef LOCAL_OPTIMIZATION_DLIB_H
#define LOCAL_OPTIMIZATION_DLIB_H

#include "optimizationProblem.h"
#include "optimizer.h"

namespace parpe {

class OptimizerDlibLineSearch : public Optimizer {
  public:
    OptimizerDlibLineSearch() = default;

    /**
     * @brief Minimize an objective function given as OptimizationProblem using
     * dlib line search algorithm
     * @param problem
     * @return Returns 0 on success.
     */

    int optimize(OptimizationProblem *problem) override;
};

} // namespace parpe

#endif
