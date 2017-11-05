#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"
#include "optimizationProblem.h"
#include "quadraticTestProblem.h"

// clang-format off
TEST_GROUP(optimizationProblem){
    void setup() {
        mock().disable();
    }

    void teardown() {
    }
};
// clang-format on

TEST(optimizationProblem, gradientChecker) {
    parpe::QuadraticTestProblem problem {};
    int numParameterIndices {1};
    int parameterIndices[numParameterIndices] {0};

    parpe::optimizationProblemGradientCheck(&problem, parameterIndices, numParameterIndices, 1e-6);
}
