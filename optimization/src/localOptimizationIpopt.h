#ifndef LOCAL_OPTIMIZATION_H
#define LOCAL_OPTIMIZATION_H

#include "optimizer.h"
#include "optimizationProblem.h"

class OptimizerIpOpt : public Optimizer {
public:
    OptimizerIpOpt();

    int optimize(OptimizationProblem *problem);
};

/**
 * @brief getLocalOptimum Get local optimum using Ipopt Optimizer
 * @param dataPath
 * @return Returns 0 on success.
 */
int getLocalOptimumIpopt(OptimizationProblem *problem);

#endif
