#ifndef QUADRATIC_TEST_PROBLEM_H
#define QUADRATIC_TEST_PROBLEM_H

#include <parpeoptimization/multiStartOptimization.h>
#include <parpeoptimization/optimizationProblem.h>

#include <gmock/gmock.h>

using ::testing::_;
using ::testing::Invoke;

namespace parpe {

/**
 * @brief The OptimizationReporterTest is a mock implementation of
 * OptimizationReporter
 */
class OptimizationReporterTest : public OptimizationReporter {
public:
    using OptimizationReporter::OptimizationReporter;

    bool starting(gsl::span<const double> parameters) const override;

    bool iterationFinished(gsl::span<const double> parameters,
                           double objectiveFunctionValue,
                           gsl::span<const double> objectiveFunctionGradient) const override;

    bool beforeCostFunctionCall(gsl::span<const double> parameters) const override;

    bool afterCostFunctionCall(gsl::span<const double> parameters,
                               double objectiveFunctionValue,
                               gsl::span<const double> objectiveFunctionGradient) const override;

    void finished(double optimalCost,
                  gsl::span<const double> parameters, int exitStatus) const override;

    bool printDebug = false;
};


class OptimizationReporterMock: public OptimizationReporter {
public:
    OptimizationReporterMock(GradientFunction *gradFun,
                         std::unique_ptr<Logger> logger)
        :OptimizationReporter(gradFun, std::make_unique<Logger>(*logger))
    {
        testRep = std::make_unique<OptimizationReporterTest>(gradFun, std::move(logger));

        ON_CALL(*this, starting(_))
            .WillByDefault(Invoke(testRep.get(), &OptimizationReporterTest::starting));
        ON_CALL(*this, iterationFinished(_, _, _))
            .WillByDefault(Invoke(testRep.get(), &OptimizationReporterTest::iterationFinished));
        ON_CALL(*this, beforeCostFunctionCall(_))
            .WillByDefault(Invoke(testRep.get(), &OptimizationReporterTest::beforeCostFunctionCall));
        ON_CALL(*this, afterCostFunctionCall(_, _, _))
            .WillByDefault(Invoke(testRep.get(), &OptimizationReporterTest::afterCostFunctionCall));
        ON_CALL(*this, finished(_, _, _))
            .WillByDefault(Invoke(testRep.get(), &OptimizationReporterTest::finished));
    }


    MOCK_CONST_METHOD1(starting, bool(gsl::span<const double> parameters));

    MOCK_CONST_METHOD3(iterationFinished,
                 bool(gsl::span<const double> parameters,
                      double objectiveFunctionValue,
                      gsl::span<const double> objectiveFunctionGradient));

    MOCK_CONST_METHOD1(beforeCostFunctionCall,
                 bool(gsl::span<const double> parameters));

    MOCK_CONST_METHOD3(afterCostFunctionCall,
                 bool(gsl::span<const double> parameters,
                      double objectiveFunctionValue,
                      gsl::span<const double> objectiveFunctionGradient));

    MOCK_CONST_METHOD3(finished,
                 void(double optimalCost, gsl::span<const double> parameters,
                      int exitStatus));

    std::unique_ptr<OptimizationReporterTest> testRep;
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
    FunctionEvaluationStatus evaluate(gsl::span<const double> parameters,
                                      double &fval,
                                      gsl::span<double> gradient,
                                      Logger *logger = nullptr,
                                      double *cpuTime = nullptr) const override;

    int numParameters() const override;
};

class QuadraticGradientFunctionMock : public GradientFunction {
public:
    QuadraticGradientFunctionMock();

    MOCK_CONST_METHOD5(evaluate_impl, FunctionEvaluationStatus(
                     gsl::span<const double> parameters,
                     double &fval,
                     gsl::span<double> gradient,
                     Logger *logger,
                     double *cpuTime));

    virtual FunctionEvaluationStatus evaluate(
            gsl::span<const double> parameters,
            double &fval,
            gsl::span<double> gradient,
            Logger *logger = nullptr,
            double *cpuTime = nullptr) const {
        return evaluate_impl(parameters, fval, gradient, logger, cpuTime);
    }

    MOCK_CONST_METHOD0(numParameters, int());


private:
    QuadraticGradientFunction fun;
};


/**
 * @brief The QuadraticTestProblem class is a test optimization problem built
 * around QuadraticGradientFunction
 */

class QuadraticTestProblem : public OptimizationProblem {
public:
    QuadraticTestProblem(std::unique_ptr<Logger> logger = std::make_unique<Logger>());
    void fillParametersMin(gsl::span<double> buffer) const override;
    void fillParametersMax(gsl::span<double> buffer) const override;

    ~QuadraticTestProblem() {
        if(!getReporterCalled && reporter)
            // manual cleanup if not passed in unique_ptr
            delete reporter;
    }
    std::unique_ptr<OptimizationReporter> getReporter() const override;

    OptimizationReporterMock *reporter;
    mutable bool getReporterCalled = false;
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
