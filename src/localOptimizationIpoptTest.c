#include "CppUTest/TestHarness_c.h"
#include "CppUTestExt/MockSupport_c.h"

#include "localOptimizationIpopt.h"
#include "testingMisc.h"
#include "tests/quadraticTestProblem.h"

extern double optimalCost;

TEST_GROUP_C_SETUP(localOptimizationIpopt) {
    optimalCost = -1e8;
}

TEST_GROUP_C_TEARDOWN(localOptimizationIpopt) {
    mock_c()->checkExpectations();
    mock_c()->clear();
}

TEST_C(localOptimizationIpopt, testOptimization) {
    OptimizationProblem *problem = getQuadraticTestProblem();
    mock_c()->expectOneCall("logFinish")->withIntParameters("exitStatus", 0);
    getLocalOptimumIpopt(problem);

    // free(problem); //?
}


