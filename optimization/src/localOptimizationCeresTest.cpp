#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"
#include "localOptimizationCeres.h"
#include "quadraticTestProblem.h"
#include "testingMisc.h"
#include <cmath>
#include <ceres/version.h>


// clang-format off
TEST_GROUP(localOptimizationCeres){
    void setup() {
        mock().clear();
    }

    void teardown() {
        mock().checkExpectations();
        mock().clear();
    }
};
// clang-format on


TEST(localOptimizationCeres, testOptimization) {
    parpe::QuadraticTestProblem problem;

    mock().expectOneCall("OptimizationReporterTest::starting");
    mock().expectOneCall("OptimizationReporterTest::finished").withIntParameter("exitStatus", 0);
    mock().ignoreOtherCalls();

    // CERES always requests gradient, this one call is for the workaround
    // below.
    //    mock().expectNCalls(1, "testObj");
    //    mock().expectNCalls(10, "testObjGrad");

    parpe::OptimizerCeres optimizer;
    auto result = optimizer.optimize(&problem);

    // check status, cost, parameter
    CHECK_EQUAL(0, std::get<0>(result));
//    DOUBLES_EQUAL(42.0, std::get<1>(result), 1e-12);
    DOUBLES_EQUAL(-1.0, std::get<2>(result).at(0), 1e-6);

    // This is a work-around for buggy ceres in ubuntu repository, which does
    // not always return the correct optimal cost
    // Needed to run on shippable.com
    auto params = std::get<2>(result);
    double optimalCost = NAN;
    problem.costFun->evaluate(params.data(), optimalCost, nullptr);
    DOUBLES_EQUAL(42.0, optimalCost, 1e-6);
}

/* Different number of calls for older versions (will fail on shippable) */
#if (CERES_VERSION_MAJOR < 1 || CERES_VERSION_MINOR < 13)
IGNORE_TEST(localOptimizationCeres, testReporterCalled) {
#else
TEST(localOptimizationCeres, testReporterCalled) {
#endif
    parpe::QuadraticTestProblem problem;
    auto o = problem.getOptimizationOptions();
    o.maxOptimizerIterations = 1;
    // to have predictable number of function calls
    o.setOption("max_num_line_search_step_size_iterations", 1);
    problem.setOptimizationOptions(o);

    // setup
    mock().expectNCalls(3, "GradientFunction::numParameters");
    mock().expectOneCall("OptimizationReporterTest::starting");

    // starting point / iteration 0
    mock().expectOneCall("OptimizationReporterTest::beforeCostFunctionCall");
    mock().expectOneCall("testObjGrad");
    mock().expectOneCall("OptimizationReporterTest::iterationFinished");
    mock().expectOneCall("OptimizationReporterTest::afterCostFunctionCall");

    // "normal" iterations
    mock().expectNCalls(o.maxOptimizerIterations, "OptimizationReporterTest::beforeCostFunctionCall");
    mock().expectNCalls(o.maxOptimizerIterations, "testObjGrad");
    mock().expectNCalls(o.maxOptimizerIterations, "OptimizationReporterTest::iterationFinished");
    mock().expectNCalls(o.maxOptimizerIterations, "OptimizationReporterTest::afterCostFunctionCall");

    mock().expectOneCall("OptimizationReporterTest::finished").ignoreOtherParameters();

    parpe::OptimizerCeres optimizer;
    optimizer.optimize(&problem);

    // don't check results. could be anywhere, due to low iteration limit
}
