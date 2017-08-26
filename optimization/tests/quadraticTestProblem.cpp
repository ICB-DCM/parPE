#include <bits/stl_tree.h>

#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"
#include "optimizationOptions.h"
#include "quadraticTestProblem.h"
#include <math.h>
#include <stdio.h>

QuadraticTestProblem::QuadraticTestProblem() {
    numOptimizationParameters = 1;
    initialParameters = new double[numOptimizationParameters]();
    parametersMin = new double[numOptimizationParameters]();
    parametersMax = new double[numOptimizationParameters]();

    optimizationOptions = new OptimizationOptions();
    optimizationOptions->maxOptimizerIterations = 12;

    optimizationOptions->optimizer = OPTIMIZER_IPOPT;

    parametersMin[0] = -1e5;
    parametersMax[0] = 1e5;
}

int QuadraticTestProblem::evaluateObjectiveFunction(const double *parameters,
                                                    double *objFunVal,
                                                    double *objFunGrad) {
    if (objFunGrad) {
        mock().actualCall("testObjGrad");

        objFunVal[0] = pow(parameters[0] + 1.0, 2) + 42.0;
        objFunGrad[0] = 2.0 * parameters[0] + 2.0;

        //        printf("g: %f %f %f\n", parameters[0], *objFunVal,
        //        objFunGrad[0]);
    } else {
        mock().actualCall("testObj");

        *objFunVal = pow(parameters[0] + 1.0, 2) + 42.0;

        //        printf("f: %f %f\n", parameters[0], *objFunVal);
    }
    return 0;
}

void QuadraticTestProblem::logOptimizerFinished(double optimalCost,
                                                const double *optimalParameters,
                                                double masterTime,
                                                int exitStatus) {
    mock().actualCall("logFinish").withIntParameter("exitStatus", exitStatus);

    this->optimalCost = optimalCost;
    this->optimalParameter = optimalParameters[0];

    //    printf("f(x) %f x %f t %f s %d\n", optimalCost, optimalParameters[0],
    //    masterTime, exitStatus);
}

QuadraticTestProblem::~QuadraticTestProblem() {
    delete[] initialParameters;
    delete[] parametersMin;
    delete[] parametersMax;
    delete optimizationOptions;
}

OptimizationProblem *
QuadraticOptimizationProblemGeneratorForMultiStart::getLocalProblemImpl(
    int multiStartIndex) {
    OptimizationProblem *problem = new QuadraticTestProblem();

    return problem;
}
