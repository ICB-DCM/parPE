#include <bits/stl_tree.h>
#include <iostream>
#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"

#include "optimizationOptions.h"
#include "testingMisc.h"
#include "localOptimizationIpopt.h"
#include "localOptimizationIpoptTNLP.h"
#include <ceres/gradient_problem_solver.h>
#include "localOptimizationCeres.h"

TEST_GROUP(optimizationOptions){
    void setup(){mock().clear();
}

void teardown() {
    mock().checkExpectations();
    mock().clear();
}
};

TEST(optimizationOptions, setGetOptionStr) {
    OptimizationOptions o;
    std::string key = "str";
    std::string expVal = "testStr";
    o.setOption(key, expVal);
    auto actVal = o.getStringOption(key);

    CHECK_EQUAL(expVal, actVal);
}

TEST(optimizationOptions, setGetOptionInt) {
    OptimizationOptions o;
    std::string key = "str";
    double expVal = 1.23;
    o.setOption(key, expVal);
    auto actVal = o.getDoubleOption(key);

    DOUBLES_EQUAL(expVal, actVal, 1e-15);
}

TEST(optimizationOptions, setGetOptionDouble) {
    OptimizationOptions o;
    std::string key = "str";
    auto expVal = 123;
    o.setOption(key, expVal);
    auto actVal = o.getIntOption(key);

    CHECK_EQUAL(expVal, actVal);
}


TEST(optimizationOptions, getNonExistingOption) {
    OptimizationOptions o;

    CHECK_THROWS(std::invalid_argument, o.getIntOption("missingKey"));
}

TEST(optimizationOptions, setIpOptOptions) {
    std::string key = "max_iter";
    int expVal = 10;

    Ipopt::SmartPtr<Ipopt::IpoptApplication> app = IpoptApplicationFactory();
    auto options = app->Options();

    OptimizationOptions o;
    o.setOption(key, expVal);
    o.for_each(setIpOptOption, &options);

    int actVal = 0;
    CHECK_EQUAL(true, options->GetIntegerValue(key, actVal, ""));
    CHECK_EQUAL(expVal, actVal);
}

TEST(optimizationOptions, setCeresOptions) {
    std::string key = "max_num_iterations";
    int expVal = 10;

    ceres::GradientProblemSolver::Options options;

    OptimizationOptions o;
    o.setOption(key, expVal);
    o.for_each(setCeresOption, &options);

    int actVal = options.max_num_iterations;

    CHECK_TRUE(options.IsValid(nullptr));
    CHECK_EQUAL(expVal, actVal);
}

