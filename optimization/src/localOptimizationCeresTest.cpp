#include <bits/stl_tree.h>

#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"

#include "localOptimizationCeres.h"
#include "testingMisc.h"
#include "quadraticTestProblem.h"
#include <math.h>


TEST_GROUP(localOptimizationCeres)
{
    void setup() {

    }

    void teardown() {
        mock().checkExpectations();
        mock().clear();
    }
};


TEST(localOptimizationCeres, testOptimization) {
    QuadraticTestProblem *problem = new QuadraticTestProblem();

    mock().expectOneCall("logFinish").withIntParameter("exitStatus", 0);
    mock().expectNCalls(0, "testObj");
    mock().expectNCalls(10, "testObjGrad");

    getLocalOptimumCeres(problem);

    DOUBLES_EQUAL(42.0, problem->optimalCost, 1e-12);
    DOUBLES_EQUAL(-1.0, problem->optimalParameter, 1e-12);

    delete problem;
}


