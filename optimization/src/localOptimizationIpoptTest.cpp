#include <bits/stl_tree.h>

#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"

#include "localOptimizationIpopt.h"
#include "optimizationOptions.h"
#include "quadraticTestProblem.h"
#include "testingMisc.h"

TEST_GROUP(localOptimizationIpopt){void setup(){}

                                   void teardown(){mock().checkExpectations();
mock().clear();
}
}
;

TEST(localOptimizationIpopt, testOptimization) {
    QuadraticTestProblem *problem = new QuadraticTestProblem();
    problem->optimizationOptions->functionTolerance = 1;

    mock().expectOneCall("logFinish").withIntParameter("exitStatus", 4);
    //    mock().expectNCalls(11, "testObj");
    //    mock().expectNCalls(12, "testObjGrad");
    mock().ignoreOtherCalls();

    OptimizerIpOpt optimizer;
    optimizer.optimize(problem);

    DOUBLES_EQUAL(42.0, problem->optimalCost, 1e-12);
    DOUBLES_EQUAL(-1.0, problem->optimalParameter, 1e-12);

    delete problem;
}
