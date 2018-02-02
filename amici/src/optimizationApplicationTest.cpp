#include <optimizationApplication.h>
#include <quadraticTestProblem.h>
#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"

// clang-format off
TEST_GROUP(optimizationApplication){
    void setup(){

    }

    void teardown(){
    }
};
// clang-format on


TEST(optimizationApplication, gradientChecker) {

        parpe::QuadraticTestProblem problem {};
        int numParameterIndices {1};
        int parameterIndices[numParameterIndices] {0};

        parpe::OptimizationApplication app;
        // TODO: currently only works for multicondition problem

        parpe::optimizationProblemGradientCheck(&problem, parameterIndices, numParameterIndices, 1e-6);
}
