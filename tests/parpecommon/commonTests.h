#include <parpecommon/parpeConfig.h>
#include <parpecommon/misc.h>
#include <parpecommon/logging.h>
#include <parpecommon/parpeException.h>
#include <parpecommon/hdf5Misc.h>

#include "testingMisc.h"

#include <cmath>
#include <vector>
#include <cstring> // strlen

#include <gtest/gtest.h>

#ifdef PARPE_ENABLE_MPI
#include <mpi.h>
#endif

using namespace parpe;

TEST(testingMisc, testTenToMinusInf) {
    ASSERT_EQ(0.0, pow(10, -INFINITY));
}

TEST(testingMisc, testWithinTolerance) {
    captureStreamToString([](){
        double atol = 0.1;
        double rtol = 0.1;

        EXPECT_TRUE(withinTolerance(1.0, 1.0, atol, rtol, 0));
        // abs false, rel true
        EXPECT_TRUE(withinTolerance(2.0, 2.15, atol, rtol, 0));
        // abs true, rel false
        EXPECT_TRUE(withinTolerance(0, 0.05, atol, rtol, 0));

        EXPECT_TRUE(withinTolerance(NAN, NAN, atol, rtol, 0));
        EXPECT_FALSE(withinTolerance(NAN, 1, atol, rtol, 0));
        EXPECT_FALSE(withinTolerance(1, NAN, atol, rtol, 0));

        EXPECT_TRUE(withinTolerance(INFINITY, INFINITY, atol, rtol, 0));
        EXPECT_FALSE(withinTolerance(1, INFINITY, atol, rtol, 0));
        EXPECT_FALSE(withinTolerance(INFINITY, 1, atol, rtol, 0));
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
        EXPECT_TRUE(r >= min && r <= max);
    }
}

TEST(commonMisc, testBacktrace) {
    std::string output = captureStreamToString([]() {
        parpe::printBacktrace(5);
    }, stderr, STDERR_FILENO);
    EXPECT_TRUE(100 < output.size());
}

TEST(commonMisc, testFilexists) {
    EXPECT_TRUE(parpe::fileExists("/"));
    EXPECT_FALSE(parpe::fileExists("/doesntExists"));
}

TEST(commonMisc, testRecursiveMkpath) {
    std::string name {"parpeTest_testRecursiveMkpath/a/b/c"};
    auto _ = gsl::finally([name] { rmdir(name.c_str()); });

    parpe::mkpathConstChar(name.c_str(), 0755);
}

TEST(commonMisc, testRandDouble) {
    const int numTests = 100;
    const double min = -1.0;
    const double max = 1.0;

    for(int i = 0; i < numTests; ++i) {
        double r = parpe::randDouble(min, max);
        EXPECT_TRUE(r >= min && r <= max);
    }
}

TEST(commonMisc, testFillArrayRandomDoubleSameInterval) {
    const int numTests = 100;
    const double min = -1.0;
    const double max = 1.0;

    double buf[numTests];
    parpe::fillArrayRandomDoubleSameInterval(min, max, buf);

    for(int i = 0; i < numTests; ++i) {
        EXPECT_TRUE(buf[i] >= min && buf[i] <= max);
    }
}

TEST(commonMisc, testFillArrayRandomDoubleIndividualInterval) {
    const int numTests = 100;
    const double min[numTests] = {-1.0, 0.0, 1.0};
    const double max[numTests] = {-0.5, 0.5, 1.5};

    double buf[numTests];
    parpe::fillArrayRandomDoubleIndividualInterval(min, max, buf);

    for(int i = 0; i < numTests; ++i) {
        EXPECT_TRUE(buf[i] >= min[i] && buf[i] <= max[i]);
    }
}


#ifdef PARPE_ENABLE_MPI
TEST(commonMisc, testMpi) {
    // Before MPI initialized
    EXPECT_EQ(-1, parpe::getMpiRank());
    EXPECT_EQ(-1, parpe::getMpiCommSize());

    // MPI initialized
    MPI_Init(0, nullptr);
    EXPECT_EQ(0, parpe::getMpiRank());
    EXPECT_EQ(1, parpe::getMpiCommSize());
    MPI_Finalize();

    // Should not make invalid calls after mpi_finalize
    EXPECT_EQ(-1, parpe::getMpiRank());
    EXPECT_EQ(-1, parpe::getMpiCommSize());
}
#endif


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
    parpe::strFormatCurrentLocaltime(gsl::make_span(buf, buflen), "abc");
    EXPECT_EQ(3UL, std::strlen(buf));
}


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

    EXPECT_TRUE(s.size() > 20);
}

TEST(logging, logProcessStats) {
    std::string s = captureStreamToString([](){
        parpe::logProcessStats();
    }, stdout);

    EXPECT_TRUE(s.size() > 200);
}


#include <parpecommon/costFunction.h>
#include <parpecommon/model.h>

TEST(costFunction, mseZero) {
    parpe::MeanSquaredError mse;
    std::vector<double> label = {1.0, 1.0};
    std::vector<double> prediction = {1.0, 1.0};
    double costExp = 0.0;
    double costAct = NAN;
    mse.evaluate(label, prediction, costAct);
    EXPECT_EQ(costExp, costAct);

    std::vector<double> predictionGradient0 = {1.0, 1.0};
    std::vector<double*> predictionGradient = { &predictionGradient0[0],
                                                &predictionGradient0[1]};

    std::vector<double> costGradientExp = {0.0};
    std::vector<double> costGradientAct = {NAN};

    mse.evaluate(label, prediction,
                 1, predictionGradient,
                 costAct, costGradientAct.data());

    EXPECT_EQ(costExp, costAct);
    EXPECT_TRUE(costGradientExp == costGradientAct);

}

TEST(costFunction, mseNonzero) {
    parpe::MeanSquaredError mse;
    std::vector<double> label = {1.0, 1.0};
    std::vector<double> prediction = {1.0, 2.0};
    double costExp = 0.5;
    double costAct = NAN;
    mse.evaluate(label, prediction, costAct);
    EXPECT_EQ(costExp, costAct);

    std::vector<double> predictionGradient0 = {5.0, 3.0};
    std::vector<double*> predictionGradient = { &predictionGradient0[0], &predictionGradient0[1]};

    std::vector<double> costGradientExp = {3.0};
    std::vector<double> costGradientAct = {NAN};

    mse.evaluate(label, prediction,
                 1, predictionGradient,
                 costAct, costGradientAct.data());

    EXPECT_EQ(costExp, costAct);
    EXPECT_TRUE(costGradientExp == costGradientAct);

}


TEST(costFunction, linearModel) {
    std::vector<double> parameters = { 3.0, 2.0 }; // x = 3.0, b = 2.0
    std::vector<std::vector<double>> features = { { 4.0 } };
    // y = A x + b = 4.0 * 3.0 + 2.0 = 14.0
    std::vector<double> outputsExp {14.0};
    std::vector<double> outputsAct(features.size(), NAN);

    LinearModel lm;
    lm.evaluate(parameters, features, outputsAct);

    EXPECT_TRUE(outputsExp == outputsAct);

    std::vector<std::vector<double>> gradExp = {{4.0, 1.0}};

    auto gradAct = std::vector<std::vector<double> >(
                features.size(),
                std::vector<double>(parameters.size(), NAN));

    lm.evaluate(parameters, features, outputsAct, gradAct);
    EXPECT_TRUE(gradExp == gradAct);
}

TEST(costFunction, linearModel2) {
    std::vector<double> parameters = { 3.0, 1.0, 2.0 };
    std::vector<std::vector<double>> features = { { 4.0, 1.0 } };
    // y = A x + b = 4.0 * 3.0 + 1.0*1.0 + 2.0 = 15.0
    std::vector<double> outputsExp {15.0};
    std::vector<double> outputsAct(features.size(), NAN);

    LinearModel lm;
    lm.evaluate(parameters, features, outputsAct);

    EXPECT_TRUE(outputsExp == outputsAct);

    std::vector<std::vector<double>> gradExp = {{4.0, 1.0, 1.0}};

    auto gradAct = std::vector<std::vector<double> >(
                features.size(),
                std::vector<double>(parameters.size(), NAN));

    lm.evaluate(parameters, features, outputsAct, gradAct);
    EXPECT_TRUE(gradExp == gradAct);
}

TEST(costFunction, linearModel3) {
    std::vector<double> parameters = { 3.0, 1.0, 2.0 };
    std::vector<std::vector<double>> features = { { 4.0, 1.0 },  { 8.0, 2.0 }};
    // y = A x + b = 4.0 * 3.0 + 1.0*1.0 + 2.0 = 15.0
    // y2= 8.0 * 3.0 + 2.0*1.0 + 2.0 = 28
    std::vector<double> outputsExp {15.0, 28.0};
    std::vector<double> outputsAct(features.size(), NAN);

    parpe::LinearModel lm;
    lm.evaluate(parameters, features, outputsAct);

    EXPECT_TRUE(outputsExp == outputsAct);

    std::vector<std::vector<double>> gradExp = {{4.0, 1.0, 1.0}, {8.0, 2.0, 1.0}};

    auto gradAct = std::vector<std::vector<double> >(
                features.size(), std::vector<double>(parameters.size(), NAN));

    lm.evaluate(parameters, features, outputsAct, gradAct);
    EXPECT_TRUE(gradExp == gradAct);
}
