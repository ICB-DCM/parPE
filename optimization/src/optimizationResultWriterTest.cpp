#include <bits/stl_tree.h>

#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"
#include "hdf5Misc.h"
#include "optimizationResultWriter.h"

TEST_GROUP(optimizationResultWriter){void setup(){parPE::initHDF5Mutex();
}

void teardown() { parPE::destroyHDF5Mutex(); }
}
;

TEST(optimizationResultWriter, testResultWriter) {
    parPE::OptimizationResultWriter w("deleteme.h5", true);

    w.setRootPath("/bla/");

    w.saveTotalCpuTime(100);

    w.flushResultWriter();

    w.logLocalOptimizerIteration(1, NULL, 0, 0, NULL, 1, 2, 3, 4, 6, 7, 8, 9,
                                 10, 11);

    w.logLocalOptimizerObjectiveFunctionEvaluation(NULL, 0, 1, NULL, 2, 3);

    w.saveLocalOptimizerResults(1, NULL, 0, 12, 0);

    w.flushResultWriter();
}
