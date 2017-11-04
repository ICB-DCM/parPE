#ifndef LOCAL_OPTIMIZATION_TOMS611_H
#define LOCAL_OPTIMIZATION_TOMS611_H

#include "optimizationProblem.h"
#include "optimizer.h"

namespace parpe {

class OptimizerToms611TrustRegionSumsl : public Optimizer {
  public:
    OptimizerToms611TrustRegionSumsl() = default;

    /**
     * @brief Minimize an objective function given as OptimizationProblem using
     * the `sumsl_` trust-region algorithm from TOMS611
     *
     * TODO: no options are specifyable for the moment
     *
     * @param problem
     * @return Returns 0 on success.
     */

    int optimize(OptimizationProblem *problem) override;
};

} // namespace parpe

#endif
