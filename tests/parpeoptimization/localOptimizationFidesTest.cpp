#include <gtest/gtest.h>

#include <parpeoptimization/localOptimizationFides.h>
#include <parpeoptimization/optimizationOptions.h>

#include "../parpecommon/testingMisc.h"
#include "quadraticTestProblem.h"

#include <fides/constants.hpp>

using ::testing::_;
using ::testing::AtLeast;
using ::testing::Eq;
using ::testing::Ne;

TEST(LocalOptimizationFides, FindsOptimum)
{
    parpe::QuadraticTestProblem problem;

    // should trigger termination
    auto xtol = 1e-8;
    // should not trigger termination
    auto fatol = 1e-16;
    auto frtol = 1e-16;
    auto gatol = 1e-16;

    auto optimization_options = problem.getOptimizationOptions();
    optimization_options.setOption("xtol", xtol);
    optimization_options.setOption("fatol", fatol);
    optimization_options.setOption("frtol", frtol);
    optimization_options.setOption("gatol", gatol);
    problem.setOptimizationOptions(optimization_options);

    EXPECT_CALL(*problem.reporter, starting(_));
    EXPECT_CALL(*problem.reporter,
                finished(_, _, static_cast<int>(fides::ExitStatus::ftol)));

    // No calls without gradient
    EXPECT_CALL(*dynamic_cast<parpe::QuadraticGradientFunctionMock*>(
                  problem.cost_fun_.get()),
                evaluate_impl(_, _, Eq(gsl::span<const double>()), _, _))
      .Times(0);
    // At least one gradient evaluation
    EXPECT_CALL(*dynamic_cast<parpe::QuadraticGradientFunctionMock*>(
                  problem.cost_fun_.get()),
                evaluate_impl(_, _, Ne(gsl::span<const double>()), _, _))
      .Times(AtLeast(1));

    parpe::OptimizerFides optimizer;
    auto [status, fval, parameters]= optimizer.optimize(&problem);

    // check status, cost, parameter
    EXPECT_EQ(0, status);
    EXPECT_NEAR(42.0, fval, 1e-12);
    EXPECT_NEAR(-1.0, parameters.at(0), xtol * 10.0);
}
