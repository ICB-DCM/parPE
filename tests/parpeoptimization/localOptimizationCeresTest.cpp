#include <gtest/gtest.h>

#include <parpeoptimization/localOptimizationCeres.h>

#include "quadraticTestProblem.h"
#include "../parpecommon/testingMisc.h"

#include <cmath>
#include <ceres/version.h>


using ::testing::_;
using ::testing::Eq;
using ::testing::Ne;
using ::testing::AtLeast;


TEST(LocalOptimizationCeres, Optimization) {
    parpe::QuadraticTestProblem problem;

    EXPECT_CALL(*problem.reporter, starting(_));
    EXPECT_CALL(*problem.reporter, finished(_, _, 0));

    // CERES always requests gradient, this one call is for the workaround
    // below.
    //    mock().expectNCalls(1, "testObj");
    //    mock().expectNCalls(10, "testObjGrad");

    parpe::OptimizerCeres optimizer;
    auto result = optimizer.optimize(&problem);

    // check status, cost, parameter
    EXPECT_EQ(0, std::get<0>(result));
    // EXPECT_NEAR(42.0, std::get<1>(result), 1e-12);
    EXPECT_NEAR(-1.0, std::get<2>(result).at(0), 1e-6);

    // This is a work-around for buggy ceres in ubuntu repository, which does
    // not always return the correct optimal cost
    // Needed to run on shippable.com
    auto params = std::get<2>(result);
    double optimalCost = NAN;
    problem.cost_fun_->evaluate(params, optimalCost, gsl::span<double>());
    EXPECT_NEAR(42.0, optimalCost, 1e-6);
}

/* Different number of calls for older versions (will fail on shippable) */
#if (CERES_VERSION_MAJOR < 1 && CERES_VERSION_MINOR < 13)
IGNORE_TEST(localOptimizationCeres, testReporterCalled) {
#else
TEST(LocalOptimizationCeres, IsReporterCalled) {
#endif
    parpe::QuadraticTestProblem problem;
    auto o = problem.getOptimizationOptions();
    o.maxOptimizerIterations = 1;
    // to have predictable number of function calls
    o.setOption("max_num_line_search_step_size_iterations", 1);
    problem.setOptimizationOptions(o);

    EXPECT_CALL(*problem.reporter, starting(_));
    EXPECT_CALL(*dynamic_cast<parpe::QuadraticGradientFunctionMock *>(problem.cost_fun_.get()),
                numParameters()).Times(3);
    EXPECT_CALL(*problem.reporter, beforeCostFunctionCall(_)).Times(1 + o.maxOptimizerIterations);
    EXPECT_CALL(*dynamic_cast<parpe::QuadraticGradientFunctionMock *>(problem.cost_fun_.get()),
                evaluate_impl(_, _, Ne(gsl::span<const double>()), _, _)).Times(1 + o.maxOptimizerIterations);
    EXPECT_CALL(*problem.reporter, iterationFinished(_, _, _)).Times(1 + o.maxOptimizerIterations);
    EXPECT_CALL(*problem.reporter, afterCostFunctionCall(_, _, _)).Times(1 + o.maxOptimizerIterations);

    EXPECT_CALL(*problem.reporter, finished(_, _, _));

    parpe::OptimizerCeres optimizer;
    optimizer.optimize(&problem);

    // don't check results. could be anywhere, due to low iteration limit
}
