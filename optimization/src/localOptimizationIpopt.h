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

void setIpOptOption(const std::pair<const std::string, const std::string> &pair, void* arg);

#endif
