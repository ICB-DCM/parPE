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
    OptimizationResultWriter w("deleteme.h5", true);

    w.rootPath = "/bla/";

    w.saveTotalCpuTime(100);

    w.flushResultWriter();

    w.logLocalOptimizerIteration(1, NULL, 0, 0, NULL, 1, 2, 3, 4, 6, 7, 8, 9, 10, 11);

    w.logLocalOptimizerObjectiveFunctionEvaluation(NULL, 0, 1, NULL, 2, 3);

    w.saveLocalOptimizerResults(1, NULL, 0, 12, 0);

    w.flushResultWriter();
}


