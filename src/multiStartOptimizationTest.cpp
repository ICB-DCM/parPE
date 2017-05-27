#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"

#include "localOptimizationIpopt.h"
#include "testingMisc.h"
#include "tests/quadraticTestProblem.h"
#include "multiStartOptimization.h"

TEST_GROUP(multiStartOptimization)
{
    void setup() {
    }

    void teardown() {

    }
};


TEST(multiStartOptimization, testMultiStartOptimization) {
    mock().disable();

    int numStarts = 1;
    runParallelMultiStartOptimization(quadraticOptimizationProblemGeneratorForMultiStart, numStarts, true);

//    mock().expectNCalls(12, "logFinish").withIntParameter("exitStatus", 0 );
//    mock().ignoreOtherCalls();
//    mock().checkExpectations();
//    mock().clear();
    mock().enable();
}


