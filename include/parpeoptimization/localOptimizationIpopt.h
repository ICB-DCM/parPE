#ifndef LOCAL_OPTIMIZATION_H
#define LOCAL_OPTIMIZATION_H

#include <parpeoptimization/optimizationProblem.h>
#include <parpeoptimization/optimizer.h>

namespace parpe {

class OptimizerIpOpt : public Optimizer {
  public:
    OptimizerIpOpt() = default;

    /**
     * @brief getLocalOptimum Get local optimum using Ipopt Optimizer
     * @param problem
     * @return Returns 0 on success.
     */

    std::tuple<int, double, std::vector<double> > optimize(OptimizationProblem *problem) override;
};

} // namespace parpe

#endif
