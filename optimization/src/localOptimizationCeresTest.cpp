#include <bits/stl_tree.h>

#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"

#include "localOptimizationCeres.h"
#include "quadraticTestProblem.h"
#include "testingMisc.h"
#include <math.h>

TEST_GROUP(localOptimizationCeres){void setup(){mock().clear();
}

void teardown() {
    mock().checkExpectations();
    mock().clear();
}
}
;

TEST(localOptimizationCeres, testOptimization) {
    QuadraticTestProblem problem;

    mock().expectOneCall("logFinish").withIntParameter("exitStatus", 0);
    // CERES always requests gradient, this one call is for the workaround
    // below.
    mock().expectNCalls(1, "testObj");
    //    mock().expectNCalls(10, "testObjGrad");
    mock().ignoreOtherCalls();

    OptimizerCeres optimizer;
    optimizer.optimize(&problem);

    // This is a work-around for buggy ceres in ubuntu repository, which does
    // not always return the correct optimal cost
    problem.evaluateObjectiveFunction(&problem.optimalParameter,
                                      &problem.optimalCost, NULL);

    DOUBLES_EQUAL(42.0, problem.optimalCost, 1e-12);
    DOUBLES_EQUAL(-1.0, problem.optimalParameter, 1e-12);
}
