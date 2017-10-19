#ifndef OPTIMIZER_H
#define OPTIMIZER_H

namespace parpe {

class OptimizationProblem;
class OptimizationOptions;

class Optimizer {
  public:
    Optimizer();

    virtual int optimize(OptimizationProblem *) = 0;

    OptimizationOptions *options = nullptr;

    virtual ~Optimizer();
};

} // namespace parpe

#endif // OPTIMIZER_H
