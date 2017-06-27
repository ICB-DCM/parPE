#include <bits/stl_tree.h>

#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"

#include "localOptimizationIpopt.h"
#include "testingMisc.h"
#include "quadraticTestProblem.h"
#include "multiStartOptimization.h"


TEST_GROUP(multiStartOptimization)
{
    void setup() {
        mock().clear();
    }

    void teardown() {
        mock().checkExpectations();
        mock().clear();
    }
};


TEST(multiStartOptimization, testMultiStartOptimization) {
    int numStarts = 4;

    // exit status may change depending on starting point -> ignore
    mock().expectNCalls(numStarts, "logFinish").ignoreOtherParameters();
    mock().ignoreOtherCalls();

    QuadraticOptimizationProblemGeneratorForMultiStart generator;
    runParallelMultiStartOptimization(&generator, numStarts, true);
}


