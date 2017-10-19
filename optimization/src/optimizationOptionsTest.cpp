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

// need prototype here, otherwise mess with headers (including ceres.h causes some errors with EIGEN)
namespace parPE {
void setCeresOption(const std::pair<const std::string, const std::string> &pair, ceres::GradientProblemSolver::Options* options);
}


TEST_GROUP(optimizationOptions){
    void setup(){mock().clear();
}

void teardown() {
    mock().checkExpectations();
    mock().clear();
}
};

TEST(optimizationOptions, setGetOptionStr) {
    parPE::OptimizationOptions o;
    std::string key = "str";
    std::string expVal = "testStr";
    o.setOption(key, expVal);
    auto actVal = o.getStringOption(key);

    CHECK_EQUAL(expVal, actVal);
}

TEST(optimizationOptions, setGetOptionInt) {
    parPE::OptimizationOptions o;
    std::string key = "str";
    double expVal = 1.23;
    o.setOption(key, expVal);
    auto actVal = o.getDoubleOption(key);

    DOUBLES_EQUAL(expVal, actVal, 1e-15);
}

TEST(optimizationOptions, setGetOptionDouble) {
    parPE::OptimizationOptions o;
    std::string key = "str";
    auto expVal = 123;
    o.setOption(key, expVal);
    auto actVal = o.getIntOption(key);

    CHECK_EQUAL(expVal, actVal);
}


TEST(optimizationOptions, getNonExistingOption) {
    parPE::OptimizationOptions o;

    CHECK_THROWS(std::invalid_argument, o.getIntOption("missingKey"));
}

TEST(optimizationOptions, setIpOptOptions) {
    std::string key = "max_iter";
    int expVal = 10;

    Ipopt::SmartPtr<Ipopt::IpoptApplication> app = IpoptApplicationFactory();
    auto options = app->Options();

    parPE::OptimizationOptions o;
    o.setOption(key, expVal);
    o.for_each<Ipopt::SmartPtr<Ipopt::OptionsList>*>(parPE::setIpOptOption, &options);

    int actVal = 0;
    CHECK_EQUAL(true, options->GetIntegerValue(key, actVal, ""));
    CHECK_EQUAL(expVal, actVal);
}

TEST(optimizationOptions, setCeresOptions) {
    std::string key = "max_num_iterations";
    int expVal = 10;

    ceres::GradientProblemSolver::Options options;

    parPE::OptimizationOptions o;
    o.setOption(key, expVal);
    o.for_each<ceres::GradientProblemSolver::Options*>(parPE::setCeresOption, &options);

    int actVal = options.max_num_iterations;

    CHECK_TRUE(options.IsValid(nullptr));
    CHECK_EQUAL(expVal, actVal);
}

