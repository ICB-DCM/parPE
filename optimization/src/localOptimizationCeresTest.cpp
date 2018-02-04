#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"
#include "localOptimizationCeres.h"
#include "quadraticTestProblem.h"
#include "testingMisc.h"
#include <cmath>


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
    DOUBLES_EQUAL(42.0, std::get<1>(result), 1e-12);
    DOUBLES_EQUAL(-1.0, std::get<2>(result).at(0), 1e-6);

    //    // This is a work-around for buggy ceres in ubuntu repository, which does
    //    // not always return the correct optimal cost
    //    problem.costFun->evaluate(&problem.optimalParameter,
    //                                      problem.optimalCost, nullptr);
    //    DOUBLES_EQUAL(42.0, problem.optimalCost, 1e-12);
    //    DOUBLES_EQUAL(-1.0, problem.optimalParameter, 1e-12);
}
