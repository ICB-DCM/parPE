#ifndef OPTIMIZER_H
#define OPTIMIZER_H

class OptimizationProblem;
class OptimizationOptions;

class Optimizer
{
public:
    Optimizer() {}

    virtual int optimize(OptimizationProblem *) = 0;

    OptimizationOptions *options = nullptr;

};

#endif // OPTIMIZER_H
