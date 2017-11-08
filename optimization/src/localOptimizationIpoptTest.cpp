#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"
#include "localOptimizationIpopt.h"
#include "optimizationOptions.h"
#include "quadraticTestProblem.h"
#include "testingMisc.h"


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
    //problem->optimizationOptions->functionTolerance = 1;

    mock().expectOneCall("logFinish").withIntParameter("exitStatus", 0);
    //    mock().expectNCalls(11, "testObj");
    //    mock().expectNCalls(12, "testObjGrad");
    mock().ignoreOtherCalls();

    parpe::OptimizerIpOpt optimizer;
    optimizer.optimize(&problem);

    DOUBLES_EQUAL(42.0, problem.optimalCost, 1e-12);
    DOUBLES_EQUAL(-1.0, problem.optimalParameter, 1e-12);
}
