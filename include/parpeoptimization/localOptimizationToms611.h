#ifndef LOCAL_OPTIMIZATION_TOMS611_H
#define LOCAL_OPTIMIZATION_TOMS611_H

#include <parpeoptimization/optimizationProblem.h>
#include <parpeoptimization/optimizer.h>

namespace parpe {

class OptimizerToms611TrustRegionSumsl : public Optimizer {
  public:
    OptimizerToms611TrustRegionSumsl() = default;

    /**
     * @brief Minimize an objective function given as OptimizationProblem using
     * the `sumsl_` trust-region algorithm from TOMS611
     * (ACM TOMS 9 (1983) 503-524, https://doi.org/10.1145/356056.356066,
     * code: http://netlib.org/toms/)
     *
     * TODO: no options are specifiable for the moment
     *
     * @param problem
     * @return Returns 0 on success.
     */

    virtual std::tuple<int, double, std::vector<double> >  optimize(OptimizationProblem *problem) override;
};

} // namespace parpe

#endif
