#include <CppUTest/TestHarness.h>
#include <CppUTestExt/MockSupport.h>
#include <localOptimizationIpopt.h>
#include <optimizationOptions.h>
#include <quadraticTestProblem.h>
#include <testingMisc.h>


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


TEST(localOptimizationIpopt, testOptimization) {
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
