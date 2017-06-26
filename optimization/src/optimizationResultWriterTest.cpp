#include <bits/stl_tree.h>

#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"
#include "optimizationResultWriter.h"
#include "hdf5Misc.h"

TEST_GROUP(optimizationResultWriter)
{
    void setup() {
        initHDF5Mutex();
    }

    void teardown() {
        destroyHDF5Mutex();
    }
};


TEST(optimizationResultWriter, testResultWriter) {
    OptimizationResultWriter w(NULL, "deleteme.h5", true);

    w.rootPath = "/bla/";

    w.saveTotalWalltime(100);

    w.flushResultWriter();

    w.logLocalOptimizerIteration(1, NULL, 2, NULL, 1, 2, 3, 4, 6, 7, 8, 9, 10, 11, 12);

    w.logLocalOptimizerObjectiveFunctionEvaluation(NULL, 1, NULL, 2, 3);

    w.logSimulation(NULL, 1, NULL, 1, 1, 2, NULL, NULL, NULL, 1, 2);

    w.saveLocalOptimizerResults(1, NULL, 12, 0);

    w.flushResultWriter();
}


