#ifndef OPTIMIZER_H
#define OPTIMIZER_H

#include <tuple>
#include <vector>

namespace parpe {

class OptimizationProblem;
class OptimizationOptions;

/**
 * @brief The Optimizer class is the interface to all batch optimizers available
 */
class Optimizer {
  public:
    virtual std::tuple<int, double, std::vector<double> > optimize(OptimizationProblem *) = 0;

    OptimizationOptions *options = nullptr;

    virtual ~Optimizer() = default;
};

} // namespace parpe

#endif // OPTIMIZER_H
