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

    parpe::OptimizationResultWriter w(tmpFilename, true, "/bla/");

    H5::H5File file(tmpFilename, H5F_ACC_RDONLY);

    w.setRootPath("/bla2");

    CHECK_TRUE(parpe::hdf5GroupExists(file.getId(), "/bla"));

    w.logOptimizerIteration(1, gsl::span<double>(), 0.0, gsl::span<double>(), 1.0);
    // should it be possible to have the same iteration twice?
    w.logOptimizerIteration(1, gsl::span<double>(), 0.0, gsl::span<double>(), 1.0);

    w.logObjectiveFunctionEvaluation(gsl::span<double>(), 1.0, gsl::span<double>(), 1, 2, 3.0);

    w.saveOptimizerResults(1.0, gsl::span<double>(), 12, 0);

    // TODO: check output

    CHECK_FALSE(remove(tmpFilename));
}
