#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"

#include <misc.h>
#include <logging.h>
#include <parpeException.h>
#include <hdf5Misc.h>
#include <mpi.h>

// clang-format off
TEST_GROUP(commonMisc){
    void setup(){

    }

    void teardown(){
        mock().checkExpectations();
        mock().clear();
    }
};
// clang-format on

TEST(commonMisc, testBacktrace) {
    parpe::printBacktrace(5);
}

TEST(commonMisc, testFilexists) {
    CHECK_TRUE(parpe::fileExists("/"));
    CHECK_FALSE(parpe::fileExists("/doesntExists"));
}

TEST(commonMisc, testCreateDirectoryIfNotExists) {
    char dir[] {"/"};
    parpe::createDirectoryIfNotExists(dir);
}

TEST(commonMisc, testRandDouble) {
    const int numTests = 100;
    const double min = -1.0;
    const double max = 1.0;

    for(int i = 0; i < numTests; ++i) {
        double r = parpe::randDouble(min, max);
        CHECK_TRUE(r >= min && r <= max);
    }
}

TEST(commonMisc, testFillArrayRandomDoubleSameInterval) {
    const int numTests = 100;
    const double min = -1.0;
    const double max = 1.0;

    double buf[numTests];
    parpe::fillArrayRandomDoubleSameInterval(min, max, numTests, buf);

    for(int i = 0; i < numTests; ++i) {
        CHECK_TRUE(buf[i] >= min && buf[i] <= max);
    }
}

TEST(commonMisc, testMpi) {
    // Before MPI initialized
    CHECK_EQUAL(-1, parpe::getMpiRank());
    CHECK_EQUAL(-1, parpe::getMpiCommSize());

    // MPI initialized
    MPI_Init(0, nullptr);
    CHECK_EQUAL(0, parpe::getMpiRank());
    CHECK_EQUAL(1, parpe::getMpiCommSize());
    MPI_Finalize();

    // Should not make invalid calls after mpi_finalize
    CHECK_EQUAL(-1, parpe::getMpiRank());
    CHECK_EQUAL(-1, parpe::getMpiCommSize());
}



