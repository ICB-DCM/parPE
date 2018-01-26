#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"
#include "hdf5Misc.h"
#include "optimizationResultWriter.h"
#include <cstdio>
#include <H5Cpp.h>

// clang-format off
TEST_GROUP(optimizationResultWriter){
    void setup() {
        parpe::initHDF5Mutex();
    }

    void teardown() {
    }
};
// clang-format on


TEST(optimizationResultWriter, testResultWriter) {
    const char* tmpFilename = "deleteme.h5";

    parpe::OptimizationResultWriter w(tmpFilename, true);

    H5::H5File file(tmpFilename, H5F_ACC_RDONLY);

    w.setRootPath("/bla/");

    CHECK_TRUE(parpe::hdf5GroupExists(file.getId(), "/bla"));

    w.logLocalOptimizerIteration(1, nullptr, 0, 0, nullptr, 1);
    // should it be possible to have the same iteration twice?
    w.logLocalOptimizerIteration(1, nullptr, 0, 0, nullptr, 1);

    w.logLocalOptimizerObjectiveFunctionEvaluation(NULL, 0, 1, NULL, 1, 2, 3);

    w.saveLocalOptimizerResults(1, NULL, 0, 12, 0);

    // TODO: check output

    CHECK_FALSE(remove(tmpFilename));
}
