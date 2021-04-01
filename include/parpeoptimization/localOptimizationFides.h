#ifndef LOCAL_OPTIMIZATION_FIDES_H
#define LOCAL_OPTIMIZATION_FIDES_H

#include <parpeoptimization/optimizationProblem.h>
#include <parpeoptimization/optimizer.h>

namespace parpe {

class OptimizerFides : public Optimizer
{
  public:
    OptimizerFides() = default;

    /**
     * @brief Minimize an objective function given as OptimizationProblem using
     * fides optimizer
     *
     * @param problem
     * @return .
     */

    std::tuple<int, double, std::vector<double>> optimize(
      OptimizationProblem* problem) override;
};

} // namespace parpe

#endif
