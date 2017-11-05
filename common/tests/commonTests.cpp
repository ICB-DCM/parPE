#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"

#include <misc.h>
#include <logging.h>
#include <parpeException.h>
#include <hdf5Misc.h>
#include <testingMisc.h>
#include <mpi.h>
#include <cmath>
#include <vector>

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

    char tmpName[TMP_MAX];
    std::tmpnam(tmpName);
    parpe::createDirectoryIfNotExists(tmpName);
    rmdir(tmpName);
}

TEST(commonMisc, testRecursiveMkpath) {
    char tmpName[TMP_MAX];
    std::tmpnam(tmpName);
    std::string name(tmpName);
    name += "/a/b/c";
    parpe::mkpathConstChar(name.c_str(), 0755);
    rmdir(name.c_str());
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

TEST(commonMisc, testFillArrayRandomDoubleIndividualInterval) {
    const int numTests = 100;
    const double min[numTests] = {-1.0, 0.0, 1.0};
    const double max[numTests] = {-0.5, 0.5, 1.5};

    double buf[numTests];
    parpe::fillArrayRandomDoubleIndividualInterval(min, max, numTests, buf);

    for(int i = 0; i < numTests; ++i) {
        CHECK_TRUE(buf[i] >= min[i] && buf[i] <= max[i]);
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


TEST(commonMisc, runInParallelAndWaitForFinish) {
    captureStreamToString([](){
        const int numThreads = 15;
        void* args[numThreads];

        parpe::runInParallelAndWaitForFinish(
                    [](void *) -> void* { return nullptr; }, args, numThreads);
    }, stdout);
}

TEST(commonMisc, strFormatCurrentLocaltime) {
    int buflen = 10;
    char buf[buflen];
    parpe::strFormatCurrentLocaltime(buf, buflen, "abc");
    CHECK_EQUAL(3, strlen(buf));
}


// clang-format off
TEST_GROUP(logging){
    void setup(){

    }

    void teardown(){
        mock().checkExpectations();
        mock().clear();
    }
};
// clang-format on


TEST(logging, printDebugInfoAndWait) {
    captureStreamToString([](){
        parpe::printDebugInfoAndWait(0);
    }, stdout);
}

TEST(logging, misc) {
    captureStreamToString([](){
        parpe::warning("bla");
        parpe::error("bla");
        parpe::logmessage(parpe::LOGLVL_ERROR, "error");
    }, stdout);
}

TEST(logging, printMPIInfo) {
    std::string s = captureStreamToString([](){
        parpe::printMPIInfo();
    }, stdout);

    CHECK_TRUE(s.size() > 20);
}

TEST(logging, logProcessStats) {
    std::string s = captureStreamToString([](){
        parpe::logProcessStats();
    }, stdout);

    CHECK_TRUE(s.size() > 200);
}


