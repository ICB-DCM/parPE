#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"
#include "optimizationProblem.h"
#include "quadraticTestProblem.h"
#include <cmath>

// clang-format off
TEST_GROUP(optimizationProblem){
    void setup() {
        mock().disable();
    }

    void teardown() {
    }
};
// clang-format on


TEST(optimizationProblem, quadraticTestFunction) {
    // Test QuadraticGradientFunction for f(-1) = 42
    parpe::QuadraticGradientFunction f {};
    double parameter = -1;
    double fValExp = 42.0;
    double gradientExp = 0;

    double fValAct = NAN;
    double gradientAct = NAN;
    f.evaluate(&parameter, fValAct, &gradientAct);

    CHECK_EQUAL(fValExp, fValAct);
    CHECK_EQUAL(gradientExp, gradientAct);
}


TEST(optimizationProblem, gradientChecker) {
    parpe::QuadraticTestProblem problem {};
    int numParameterIndices {1};
    int parameterIndices[numParameterIndices] {0};

    parpe::optimizationProblemGradientCheck(&problem, parameterIndices, numParameterIndices, 1e-6);
}
