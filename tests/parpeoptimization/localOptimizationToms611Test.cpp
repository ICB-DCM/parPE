#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"

#include "optimizationOptions.h"
#include "testingMisc.h"
#include <cmath>
#include <iostream>
#include <localOptimizationToms611.h>
#include <quadraticTestProblem.h>

#include <toms611.h> // must include last due to some interference with other headers

// clang-format off
TEST_GROUP(localOptimizationToms611){
    void setup() {
        mock().clear();
    }

    void teardown() {
        mock().checkExpectations();
        mock().clear();
    }
};
// clang-format on



void calcf(integer &n, doublereal *x, integer &nf, doublereal &f,
           integer *uiparm, doublereal *urparm, void *ufparm) {
    f = pow(x[0] + 1.0, 2) + 42.0;
    std::cout<<nf<<": "<<f<<std::endl;

}

void calcg(integer &n, doublereal *x, integer &nf, doublereal *g,
           integer *uiparm, doublereal *urparm, void *ufparm) {
    g[0] = 2.0 * x[0] + 2.0;
    std::cout<<nf<<"g: "<<g[0]<<std::endl;
}

TEST(localOptimizationToms611, testOptimization) {
    constexpr integer numOptimizationVariables = 1;

    integer liv = toms611_sumsl_iv_min_length;
    integer iv[liv];
    iv[0] = 0; // fresh start, make sumsl_ call deflt_

    integer lv = toms611_sumsl_v_min_length(numOptimizationVariables);
    doublereal v[lv];

    doublereal scaling[numOptimizationVariables] = {1};
    doublereal startingPoint[numOptimizationVariables] = {1234567890};

    integer numOptimizationVariables2 = numOptimizationVariables;
    sumsl_(numOptimizationVariables2,
           scaling,
           startingPoint,
           reinterpret_cast<S_fp>(calcf), (S_fp)calcg,
           iv, liv,
           lv, v,
           nullptr, nullptr, nullptr);

    CHECK_EQUAL(relative_function_convergence, iv[0]);
    DOUBLES_EQUAL(-1.0, startingPoint[0], 1e-7);

    //std::cout<<std::endl;
    //std::cout<<"iv[0] = "<<iv[0]<<std::endl;
    //std::cout<<"Final x = "<<startingPoint[0]<<std::endl;
}


TEST(localOptimizationToms611, testOptimizationGetlocalOptimum) {
    parpe::QuadraticTestProblem problem;

    mock().expectOneCall("OptimizationReporterTest::starting");
    mock().expectOneCall("OptimizationReporterTest::finished").withIntParameter("exitStatus", x_and_relative_function_convergence);
//    mock().expectNCalls(2, "testObj");
//    mock().expectNCalls(2, "testObjGrad");
    mock().ignoreOtherCalls();

    parpe::OptimizerToms611TrustRegionSumsl optimizer;
    auto result = optimizer.optimize(&problem);

    // check status, cost, parameter
    CHECK_EQUAL(0, std::get<0>(result));
    DOUBLES_EQUAL(42.0, std::get<1>(result), 1e-12);
    DOUBLES_EQUAL(-1.0, std::get<2>(result).at(0), 1e-8); // TODO adapt to optimizer tolerances
}


TEST(localOptimizationToms611, testReporterCalled) {
    const int allowedLineSearchSteps = 3; // only true when iteration limit == 1
    parpe::QuadraticTestProblem problem;
    auto o = problem.getOptimizationOptions();
    o.maxOptimizerIterations = 1;
    o.setOption("mxfcal", allowedLineSearchSteps + 1); // +1 for initial function evaluation
    problem.setOptimizationOptions(o);

    // iteration 0
    mock().expectNCalls(3, "GradientFunction::numParameters");
    mock().expectOneCall("OptimizationReporterTest::starting");
    mock().expectNCalls(2, "OptimizationReporterTest::beforeCostFunctionCall");
    // one should be enough:
    mock().expectNCalls(1, "testObj");
    mock().expectNCalls(1, "testObjGrad");
    mock().expectNCalls(2, "OptimizationReporterTest::afterCostFunctionCall");
    mock().expectOneCall("OptimizationReporterTest::iterationFinished");

    // others
    mock().expectNCalls(o.maxOptimizerIterations * allowedLineSearchSteps + 1, "OptimizationReporterTest::beforeCostFunctionCall");
    mock().expectNCalls(o.maxOptimizerIterations * allowedLineSearchSteps, "testObj");
    mock().expectNCalls(o.maxOptimizerIterations, "testObjGrad");
    mock().expectNCalls(o.maxOptimizerIterations, "OptimizationReporterTest::iterationFinished");
    mock().expectNCalls(o.maxOptimizerIterations * allowedLineSearchSteps + 1, "OptimizationReporterTest::afterCostFunctionCall");

    mock().expectOneCall("OptimizationReporterTest::finished").ignoreOtherParameters();

    parpe::OptimizerToms611TrustRegionSumsl optimizer;
    optimizer.optimize(&problem);

    // don't check results. could be anywhere, due to low iteration limit
}
