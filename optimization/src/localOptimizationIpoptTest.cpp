#include <bits/stl_tree.h>

#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"

#include "localOptimizationIpopt.h"
#include "optimizationOptions.h"
#include "quadraticTestProblem.h"
#include "testingMisc.h"

TEST_GROUP(localOptimizationIpopt){void setup(){}

                                   void teardown(){mock().checkExpectations();
mock().clear();
}
}
;

extern "C" {
  void deflt_(int &alg, int *iv, int &liv, int &lv, double *v);

  void sumsl_(
    int &n, double *d, double *x,
    void (*calcf)(int &, double *, int &, double &, int *, double *, void *),
    void (*calcg)(int &, double *, int &, double *, int *, double *, void *),
    int *iv, int &liv, int &lv, double *v,
    int *uiparm, double *urparm, void *ufparm);
}

TEST(localOptimizationIpopt, testOptimization) {
    parpe::QuadraticTestProblem problem;
    //problem->optimizationOptions->functionTolerance = 1;

    mock().expectOneCall("logFinish").withIntParameter("exitStatus", 0);
    //    mock().expectNCalls(11, "testObj");
    //    mock().expectNCalls(12, "testObjGrad");
    mock().ignoreOtherCalls();

    parpe::OptimizerIpOpt optimizer;
    optimizer.optimize(&problem);

    DOUBLES_EQUAL(42.0, problem.optimalCost, 1e-12);
    DOUBLES_EQUAL(-1.0, problem.optimalParameter, 1e-12);
}
