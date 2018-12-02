#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"

#include "parpeConfig.h"
#include "multiStartOptimization.h"
#include "quadraticTestProblem.h"
#include "testingMisc.h"

#ifdef PARPE_ENABLE_IPOPT
#include "localOptimizationIpopt.h"

// clang-format off
TEST_GROUP(multiStartOptimization){
    void setup(){
        // Disable mock for multi-threaded test; mocking does not seem to be thread-safe
        mock().disable();
    }

    void teardown(){
    }
};
// clang-format on


TEST(multiStartOptimization, testMultiStartOptimization) {
    int numStarts = 10;

    // exit status may change depending on starting point -> ignore
    // mock().expectNCalls(numStarts, "logFinish").ignoreOtherParameters();
    // mock().ignoreOtherCalls();

    parpe::QuadraticOptimizationMultiStartProblem ms(numStarts, true);
    parpe::MultiStartOptimization optimizer(ms);

    optimizer.runMultiThreaded();

    // TODO: check calls
    optimizer.runSingleThreaded();
}

// TODO: test retry on error
#endif
