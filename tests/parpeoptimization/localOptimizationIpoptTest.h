#include <gtest/gtest.h>

#include <parpeoptimization/localOptimizationIpopt.h>
#include <parpeoptimization/optimizationOptions.h>

#include "quadraticTestProblem.h"
#include "../parpecommon/testingMisc.h"

using ::testing::_;
using ::testing::Eq;
using ::testing::Ne;
using ::testing::AtLeast;

TEST(localOptimizationIpopt, testOptimizationResult) {
    parpe::QuadraticTestProblem problem;

    EXPECT_CALL(*problem.reporter, starting(_));
    EXPECT_CALL(*problem.reporter, finished(_, _, 0));

    EXPECT_CALL(*dynamic_cast<parpe::QuadraticGradientFunctionMock *>(problem.cost_fun_.get()),
                evaluate_impl(_, _, Eq(gsl::span<const double>()), _, _)).Times(AtLeast(1));
    EXPECT_CALL(*dynamic_cast<parpe::QuadraticGradientFunctionMock *>(problem.cost_fun_.get()),
                evaluate_impl(_, _, Ne(gsl::span<const double>()), _, _)).Times(AtLeast(1));

    // TODO mock().ignoreOtherCalls();

    parpe::OptimizerIpOpt optimizer;
    auto result = optimizer.optimize(&problem);

    // check status, cost, parameter
    EXPECT_EQ(0, std::get<0>(result));
    EXPECT_NEAR(42.0, std::get<1>(result), 1e-12);
    EXPECT_NEAR(-1.0, std::get<2>(result).at(0), 1e-12);
}

TEST(localOptimizationIpopt, testReporterCalled) {
    parpe::QuadraticTestProblem problem;
    auto o = problem.getOptimizationOptions();
    o.maxOptimizerIterations = 1;
    // to have predictable number of function calls
    o.setOption("accept_every_trial_step", "yes");
    problem.setOptimizationOptions(o);

    EXPECT_CALL(*problem.reporter, starting(_));

    EXPECT_CALL(*problem.reporter, beforeCostFunctionCall(_)).Times(3 + o.maxOptimizerIterations * 2);
    EXPECT_CALL(*dynamic_cast<parpe::QuadraticGradientFunctionMock *>(problem.cost_fun_.get()),
                evaluate_impl(_, _, Ne(gsl::span<const double>()), _, _)).Times(1 + o.maxOptimizerIterations);
    EXPECT_CALL(*dynamic_cast<parpe::QuadraticGradientFunctionMock *>(problem.cost_fun_.get()),
                evaluate_impl(_, _, Eq(gsl::span<const double>()), _, _)).Times(o.maxOptimizerIterations);
    EXPECT_CALL(*dynamic_cast<parpe::QuadraticGradientFunctionMock *>(problem.cost_fun_.get()),
                numParameters()).Times(2 + 0*o.maxOptimizerIterations);
    EXPECT_CALL(*problem.reporter, iterationFinished(_, _, _)).Times(1 + o.maxOptimizerIterations);
    EXPECT_CALL(*problem.reporter, afterCostFunctionCall(_, _, _)).Times(3 + o.maxOptimizerIterations * 2);

    EXPECT_CALL(*problem.reporter, finished(_, _, _));

    parpe::OptimizerIpOpt optimizer;
    optimizer.optimize(&problem);

    // don't check results. could be anywhere, due to low iteration limit
}

