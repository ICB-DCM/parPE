#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"

#include "localOptimizationCeres.h"
#include "testingMisc.h"
#include "tests/quadraticTestProblem.h"
#include <math.h>

extern double quadraticTestProblemOptimalCost, quadraticTestProblemOptimalParameter;


TEST_GROUP(localOptimizationCeres)
{
    void setup() {
        quadraticTestProblemOptimalCost = NAN;
        quadraticTestProblemOptimalParameter = NAN;

    }

    void teardown() {
        mock().checkExpectations();
        mock().clear();

    }
};


TEST(localOptimizationCeres, testOptimization) {
    OptimizationProblem *problem = getQuadraticTestProblem();

    mock().expectOneCall("logFinish").withIntParameter("exitStatus", 0);
    mock().expectNCalls(0, "testObj");
    mock().expectNCalls(10, "testObjGrad");

    getLocalOptimumCeres(problem);
    DOUBLES_EQUAL(42.0, quadraticTestProblemOptimalCost, 1e-12);
    DOUBLES_EQUAL(-1.0, quadraticTestProblemOptimalParameter, 1e-12);

    delete problem;
}


