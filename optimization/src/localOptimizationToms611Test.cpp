#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"

#include "optimizationOptions.h"
#include "testingMisc.h"
#include <cmath>
#include <iostream>
#include <localOptimizationToms611.h>
#include <quadraticTestProblem.h>

#include <toms611.h>

// clang-format off
TEST_GROUP(localOptimizationToms611){
    void setup(){
    }

    void teardown(){
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

IGNORE_TEST(localOptimizationToms611, testOptimization) {
    integer numOptimizationVariables = 1;

    integer liv = toms611_sumsl_iv_min_length;
    integer iv[liv];
    iv[0] = 0; // fresh start, make sumsl_ call deflt_

    integer lv = toms611_sumsl_v_min_length(numOptimizationVariables);
    doublereal v[lv];

    doublereal scaling[numOptimizationVariables] = {1};
    doublereal startingPoint[numOptimizationVariables] = {1234567890};

    sumsl_(numOptimizationVariables,
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
    problem.printDebug = true;

    mock().expectOneCall("logFinish").withIntParameter("exitStatus", x_and_relative_function_convergence);
    mock().expectNCalls(2, "testObj");
    mock().expectNCalls(2, "testObjGrad");
    mock().ignoreOtherCalls();

    parpe::OptimizerToms611TrustRegionSumsl optimizer;
    optimizer.optimize(&problem);

    DOUBLES_EQUAL(42.0, problem.optimalCost, 1e-12);
    DOUBLES_EQUAL(-1.0, problem.optimalParameter, 1e-12);
}
