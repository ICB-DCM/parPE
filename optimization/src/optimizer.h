#ifndef OPTIMIZER_H
#define OPTIMIZER_H

namespace parPE {

class OptimizationProblem;
class OptimizationOptions;

class Optimizer {
  public:
    Optimizer();

    virtual int optimize(OptimizationProblem *) = 0;

    OptimizationOptions *options = nullptr;

    virtual ~Optimizer();
};

} // namespace parPE

#endif // OPTIMIZER_H
