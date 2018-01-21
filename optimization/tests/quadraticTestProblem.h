#ifndef QUADRATIC_TEST_PROBLEM_H
#define QUADRATIC_TEST_PROBLEM_H

#include "multiStartOptimization.h"
#include "optimizationProblem.h"

namespace parpe {

/**
 * @brief The QuadraticGradientFunction class is a simple function for testing
 * the optimization framework.
 *
 * Represents the function
 *   f(x) = (x + 1)^2 + 42 = x^2 + 2x + 43
 *   for x in R
 * with gradient
 *   f'(x) = 2x + 2
 * and minimum 42 at x = -1
 */

class QuadraticGradientFunction : public GradientFunction {
public:
    FunctionEvaluationStatus evaluate(
            const double* const parameters,
            double &fval,
            double* gradient) const override;

    int numParameters() const override;
};

class QuadraticTestProblem : public OptimizationProblem {
  public:
    QuadraticTestProblem();
    /*
     * Test Problem for minimization
     * f (x) = (x+1)^2 + 42 = x^2  + 2x + 1 + 42
     * f'(x) = 2x + 2 = 0 <=> x = -1
     * f(-1) = 42
     */

    int evaluateObjectiveFunction(const double *parameters, double *objFunVal,
                                  double *objFunGrad) override;

    void logOptimizerFinished(double optimalCost,
                              const double *optimalParameters,
                              double masterTime, int exitStatus) override;

    double optimalCost;
    double optimalParameter;
    bool printDebug = false;
};

class QuadraticOptimizationMultiStartProblem : public MultiStartOptimization {
  public:
    using MultiStartOptimization::MultiStartOptimization;
    std::unique_ptr<OptimizationProblem> getLocalProblemImpl(int multiStartIndex) override;
};

} // namespace parpe

#endif
