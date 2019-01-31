#include <CppUTest/TestHarness.h>
#include <CppUTestExt/MockSupport.h>

#include <parpeoptimization/localOptimizationIpopt.h>
#include <parpeoptimization/optimizationOptions.h>

#include "quadraticTestProblem.h"
#include "../parpecommon/testingMisc.h"


// clang-format off
TEST_GROUP(localOptimizationIpopt){
    void setup() {
        mock().clear();
    }

    void teardown() {
        mock().checkExpectations();
        mock().clear();
    }
};
// clang-format on


TEST(localOptimizationIpopt, testOptimizationResult) {
    parpe::QuadraticTestProblem problem;

    // mock().expectOneCall("OptimizationReporterTest::starting");
    mock().expectOneCall("OptimizationReporterTest::finished").withIntParameter("exitStatus", 0);
    //    mock().expectNCalls(11, "testObj");
    //    mock().expectNCalls(12, "testObjGrad");
    mock().ignoreOtherCalls();

    parpe::OptimizerIpOpt optimizer;
    auto result = optimizer.optimize(&problem);

    // check status, cost, parameter
    CHECK_EQUAL(0, std::get<0>(result));
    DOUBLES_EQUAL(42.0, std::get<1>(result), 1e-12);
    DOUBLES_EQUAL(-1.0, std::get<2>(result).at(0), 1e-12);
}

TEST(localOptimizationIpopt, testReporterCalled) {
    parpe::QuadraticTestProblem problem;
    auto o = problem.getOptimizationOptions();
    o.maxOptimizerIterations = 1;
    // to have predictable number of function calls
    o.setOption("accept_every_trial_step", "yes");
    problem.setOptimizationOptions(o);

    // iteration 0
    mock().expectNCalls(3, "GradientFunction::numParameters");
    mock().expectOneCall("OptimizationReporterTest::starting");
    mock().expectNCalls(3, "OptimizationReporterTest::beforeCostFunctionCall");
    // one should be enough:
    mock().expectNCalls(1, "testObjGrad");
    mock().expectNCalls(3, "OptimizationReporterTest::afterCostFunctionCall");
    mock().expectOneCall("OptimizationReporterTest::iterationFinished");

    // others
    mock().expectNCalls(o.maxOptimizerIterations * 2, "OptimizationReporterTest::beforeCostFunctionCall");
    mock().expectNCalls(o.maxOptimizerIterations, "testObj");
    mock().expectNCalls(o.maxOptimizerIterations, "testObjGrad");
    mock().expectNCalls(0*o.maxOptimizerIterations, "GradientFunction::numParameters");
    mock().expectNCalls(o.maxOptimizerIterations, "OptimizationReporterTest::iterationFinished");
    mock().expectNCalls(o.maxOptimizerIterations * 2, "OptimizationReporterTest::afterCostFunctionCall");

    mock().expectOneCall("OptimizationReporterTest::finished").ignoreOtherParameters();

    parpe::OptimizerIpOpt optimizer;
    optimizer.optimize(&problem);

    // don't check results. could be anywhere, due to low iteration limit
}

