#ifndef QUADRATIC_TEST_PROBLEM_H
#define QUADRATIC_TEST_PROBLEM_H

#include "multiStartOptimization.h"
#include "optimizationProblem.h"

namespace parpe {

/**
 * @brief The OptimizationReporterTest is a mock implementation of OptimizationReporter
 */

class OptimizationReporterTest : public OptimizationReporter {
    using OptimizationReporter::OptimizationReporter;

    virtual bool starting(gsl::span<const double> parameters) const override;

    virtual bool iterationFinished(gsl::span<const double> parameters,
                                   double objectiveFunctionValue,
                                   gsl::span<const double> objectiveFunctionGradient) const override;

    virtual bool beforeCostFunctionCall(gsl::span<const double> parameters) const override;

    virtual bool afterCostFunctionCall(gsl::span<const double> parameters,
                                       double objectiveFunctionValue,
                                       gsl::span<const double> objectiveFunctionGradient) const override;

    virtual void finished(double optimalCost,
                          gsl::span<const double> parameters, int exitStatus) const override;

    bool printDebug = false;
};



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
 *
 * f (x) = (x+1)^2 + 42 = x^2  + 2x + 1 + 42
 * f'(x) = 2x + 2 = 0 <=> x = -1
 * f(-1) = 42
 */

class QuadraticGradientFunction : public GradientFunction {
public:
    FunctionEvaluationStatus evaluate(
            gsl::span<const double> parameters,
            double &fval,
            gsl::span<double> gradient) const override;

    int numParameters() const override;
};



/**
 * @brief The QuadraticTestProblem class is a test optimization problem built around QuadraticGradientFunction
 */

class QuadraticTestProblem : public OptimizationProblem {
public:
    QuadraticTestProblem();
    void fillParametersMin(gsl::span<double> buffer) const override;
    void fillParametersMax(gsl::span<double> buffer) const override;

    std::unique_ptr<OptimizationReporter> getReporter() const override;
};


class QuadraticOptimizationMultiStartProblem : public MultiStartOptimizationProblem {
public:
    QuadraticOptimizationMultiStartProblem(int numberOfStarts, bool restartOnFailure = false)
        : numberOfStarts(numberOfStarts), restartOnFailure_(restartOnFailure)
    {
        QuadraticTestProblem p;
        options = p.getOptimizationOptions();
    }

    std::unique_ptr<OptimizationProblem> getLocalProblem(int multiStartIndex) const override;

    int getNumberOfStarts() const override { return numberOfStarts; }

    virtual bool restartOnFailure() const override { return restartOnFailure_; }

    OptimizationOptions options;

private:

    int numberOfStarts = 1;
    bool restartOnFailure_ = false;
};

} // namespace parpe

#endif
