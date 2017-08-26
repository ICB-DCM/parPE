#ifndef LOCAL_OPTIMIZATION_H
#define LOCAL_OPTIMIZATION_H

#include "optimizationProblem.h"
#include "optimizer.h"

class OptimizerIpOpt : public Optimizer {
  public:
    OptimizerIpOpt();

    /**
     * @brief getLocalOptimum Get local optimum using Ipopt Optimizer
     * @param problem
     * @return Returns 0 on success.
     */

    int optimize(OptimizationProblem *problem) override;
};

#endif
