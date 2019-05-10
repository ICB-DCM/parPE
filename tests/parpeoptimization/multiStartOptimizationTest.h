#include <gtest/gtest.h>

#include <parpecommon/parpeConfig.h>
#include <parpeoptimization/multiStartOptimization.h>

#include "quadraticTestProblem.h"
#include "../parpecommon/testingMisc.h"

#ifdef PARPE_ENABLE_IPOPT
#include <parpeoptimization/localOptimizationIpopt.h>

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
