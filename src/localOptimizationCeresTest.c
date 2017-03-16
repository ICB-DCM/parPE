#include "CppUTest/TestHarness_c.h"
#include "CppUTestExt/MockSupport_c.h"

#include "localOptimizationCeres.hpp"
#include "testingMisc.h"
#include "tests/quadraticTestProblem.h"

extern double quadraticTestProblemOptimalCost, quadraticTestProblemOptimalParameter;

TEST_GROUP_C_SETUP(localOptimizationCeres) {
    quadraticTestProblemOptimalCost = -1e8;
    quadraticTestProblemOptimalParameter = -1e8;
}

TEST_GROUP_C_TEARDOWN(localOptimizationCeres) {
    mock_c()->checkExpectations();
    mock_c()->clear();
}


TEST_C(localOptimizationCeres, testOptimization) {
    OptimizationProblem *problem = getQuadraticTestProblem();

    mock_c()->expectOneCall("logFinish")->withIntParameters("exitStatus", 0);
    mock_c()->expectNCalls(0, "testObj");
    mock_c()->expectNCalls(10, "testObjGrad");

    getLocalOptimumCeres(problem);
    // TODO: disabled due to ceres bug?!
    // CHECK_EQUAL_C_REAL(42.0, quadraticTestProblemOptimalCost, 1e-12);
    CHECK_EQUAL_C_REAL(-1.0, quadraticTestProblemOptimalParameter, 1e-12);
    // free(problem); //?
}


