#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"

#include "optimizationOptions.h"
#include "testingMisc.h"
#include <cmath>
#include <iostream>
#include <quadraticTestProblem.h>

#include <localOptimizationFsqp.h>

// clang-format off
TEST_GROUP(localOptimizationFsqp){
    void setup() {
        mock().clear();
    }

    void teardown() {
        mock().checkExpectations();
        mock().clear();
    }
};
// clang-format on


TEST(localOptimizationFsqp, testOptimizationGetlocalOptimum) {
    parpe::QuadraticTestProblem problem;

    mock().ignoreOtherCalls();
    parpe::OptimizerFsqp optimizer;
    //auto result = optimizer.optimize(&problem);
    auto result = optimizer.optimize(&problem);

    // check status, cost, parameter
//    CHECK_EQUAL(0, std::get<0>(result));
    DOUBLES_EQUAL(42.0, std::get<1>(result), 1e-12);
    DOUBLES_EQUAL(-1.0, std::get<2>(result).at(0), 1e-8); // TODO adapt to optimizer tolerances
}



TEST(localOptimizationFsqp, testParallelMultistart) {
    /* Test if thread-safe
     * Test with:
     *    while ./build/optimization/tests/unittests_optimization_fsqp; do :; done
     */

    mock().disable(); // mock() is not thread-safe

    constexpr int numStarts {10};
    parpe::QuadraticOptimizationMultiStartProblem msp(numStarts, false);
    msp.options.optimizer = parpe::optimizerName::OPTIMIZER_FSQP;

    parpe::MultiStartOptimization mso(msp);
    mso.runMultiThreaded();

    mock().enable();
}



TEST(localOptimizationFsqp, testReporterCalled) {
    parpe::QuadraticTestProblem problem;
    auto o = problem.getOptimizationOptions();
    o.maxOptimizerIterations = 2;

    problem.setOptimizationOptions(o);

    // setup
    mock().expectNCalls(3, "GradientFunction::numParameters");
    mock().expectOneCall("OptimizationReporterTest::starting");

    // starting point / iteration 0
    mock().expectNCalls(2, "OptimizationReporterTest::beforeCostFunctionCall");
    mock().expectNCalls(1, "testObjGrad");
    mock().expectOneCall("OptimizationReporterTest::iterationFinished");
    mock().expectNCalls(2, "OptimizationReporterTest::afterCostFunctionCall");

    // "normal" iterations
    // "before" and "after" are called for f and g, but g is already cached
    mock().expectNCalls(o.maxOptimizerIterations * 2, "OptimizationReporterTest::beforeCostFunctionCall");
    mock().expectNCalls(o.maxOptimizerIterations, "testObjGrad");
    mock().expectNCalls(o.maxOptimizerIterations * 2, "OptimizationReporterTest::afterCostFunctionCall");
    mock().expectNCalls(o.maxOptimizerIterations, "OptimizationReporterTest::iterationFinished");

    mock().expectOneCall("OptimizationReporterTest::finished").ignoreOtherParameters();

    parpe::OptimizerFsqp optimizer;
    optimizer.optimize(&problem);

    // don't check results. could be anywhere, due to low iteration limit
}

