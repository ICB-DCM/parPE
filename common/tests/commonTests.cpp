#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"

#include <misc.h>
#include <logging.h>
#include <parpeException.h>
#include <hdf5Misc.h>
#include <testingMisc.h>
#include <mpi.h>
#include <cmath>

// clang-format off
TEST_GROUP(testingMisc){
    void setup(){

    }

    void teardown(){
        mock().checkExpectations();
        mock().clear();
    }
};
// clang-format on

TEST(testingMisc, testWithinTolerance) {
    captureStreamToString([](){
        double atol = 0.1;
        double rtol = 0.1;

        CHECK_TRUE(withinTolerance(1.0, 1.0, atol, rtol, 0));
        CHECK_TRUE(withinTolerance(2.0, 2.15, atol, rtol, 0)); // abs false, rel true
        CHECK_TRUE(withinTolerance(0, 0.05, atol, rtol, 0)); // abs true, rel false


        CHECK_TRUE(withinTolerance(NAN, NAN, atol, rtol, 0));
        CHECK_FALSE(withinTolerance(NAN, 1, atol, rtol, 0));
        CHECK_FALSE(withinTolerance(1, NAN, atol, rtol, 0));

        CHECK_TRUE(withinTolerance(INFINITY, INFINITY, atol, rtol, 0));
        CHECK_FALSE(withinTolerance(1, INFINITY, atol, rtol, 0));
        CHECK_FALSE(withinTolerance(INFINITY, 1, atol, rtol, 0));
    }, stderr, STDERR_FILENO);
}

TEST(testingMisc, testCheckEqualArray) {
    const double expected[] = {1.0, 2.0, 3.0};
    const double actual[] = {1.0, 2.0, 3.0};

    checkEqualArray(expected, actual, 3, 1e-16, 1e-16);
    checkEqualArray(nullptr, nullptr, 3, 1e-16, 1e-16);
}

TEST(testingMisc, testRandInt) {
    const int numTests = 100;
    const int min = -1;
    const int max = 1;

    for(int i = 0; i < numTests; ++i) {
        int r = randInt(min, max);
        CHECK_TRUE(r >= min && r <= max);
    }
}


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
    std::string output = captureStreamToString([]() {
        parpe::printBacktrace(5);
    }, stderr, STDERR_FILENO);
    CHECK_TRUE(100 < output.size());
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



