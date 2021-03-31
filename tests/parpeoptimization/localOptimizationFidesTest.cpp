#include <gtest/gtest.h>

#include <parpeoptimization/localOptimizationFides.h>
#include <parpeoptimization/optimizationOptions.h>

#include "quadraticTestProblem.h"
#include "../parpecommon/testingMisc.h"

using ::testing::_;
using ::testing::Eq;
using ::testing::Ne;
using ::testing::AtLeast;

TEST(localOptimizationFides, testOptimizationResult) {
    parpe::QuadraticTestProblem problem;

    EXPECT_CALL(*problem.reporter, starting(_));
    EXPECT_CALL(*problem.reporter, finished(_, _, 0));

    EXPECT_CALL(*dynamic_cast<parpe::QuadraticGradientFunctionMock *>(problem.cost_fun_.get()),
                evaluate_impl(_, _, Eq(gsl::span<const double>()), _, _)).Times(AtLeast(1));
    EXPECT_CALL(*dynamic_cast<parpe::QuadraticGradientFunctionMock *>(problem.cost_fun_.get()),
                evaluate_impl(_, _, Ne(gsl::span<const double>()), _, _)).Times(AtLeast(1));

    // TODO mock().ignoreOtherCalls();

    parpe::OptimizerFides optimizer;
    auto result = optimizer.optimize(&problem);

    // check status, cost, parameter
    EXPECT_EQ(0, std::get<0>(result));
    EXPECT_NEAR(42.0, std::get<1>(result), 1e-12);
    EXPECT_NEAR(-1.0, std::get<2>(result).at(0), 1e-12);
}
