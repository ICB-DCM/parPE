#include <gtest/gtest.h>

#include <parpecommon/hdf5Misc.h>
#include <parpeoptimization/optimizationResultWriter.h>

#include <cstdio>
#include <H5Cpp.h>

TEST(optimizationResultWriter, testResultWriter) {
    const char* tmpFilename = "deleteme.h5";

    parpe::OptimizationResultWriter w(tmpFilename, true, "/bla/");

    H5::H5File file(tmpFilename, H5F_ACC_RDONLY);

    EXPECT_TRUE(parpe::hdf5GroupExists(file.getId(), "/bla"));

    w.setRootPath("/bla2");

    w.logOptimizerIteration(1, gsl::span<double>(), 0.0, gsl::span<double>(),
                            1.0, 2.0);

    // should it be possible to have the same iteration twice?
    w.logOptimizerIteration(1, gsl::span<double>(), 0.0, gsl::span<double>(),
                            1.0, 2.0);

    w.logObjectiveFunctionEvaluation(gsl::span<double>(), 1.0,
                                     gsl::span<double>(), 1, 2, 3.0);

    w.saveOptimizerResults(1.0, gsl::span<double>(), 12.0, 17.0, 0);

    // TODO: check output

    EXPECT_FALSE(remove(tmpFilename));
}
