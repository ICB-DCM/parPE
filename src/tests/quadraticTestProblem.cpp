#include "CppUTest/TestHarness_c.h"
#include "CppUTestExt/MockSupport_c.h"

#include "quadraticTestProblem.h"
#include <math.h>
#include <stdio.h>

QuadraticTestProblem::QuadraticTestProblem()
{
    numOptimizationParameters = 1;
    initialParameters = new double[numOptimizationParameters];
    parametersMin = new double[numOptimizationParameters];
    parametersMax = new double[numOptimizationParameters];
    maxOptimizerIterations = 12;


    initialParameters[0] = -100;
    parametersMin[0] = -1e5;
    parametersMax[0] = 1e5;
}

int QuadraticTestProblem::evaluateObjectiveFunction(const double *parameters, double *objFunVal, double *objFunGrad)
{
    if(objFunGrad) {
        mock_c()->actualCall("testObjGrad");

        objFunVal[0] = pow(parameters[0] + 1.0, 2) + 42.0;
        objFunGrad[0] = 2.0 * parameters[0] + 2.0;

//        printf("g: %f %f %f\n", parameters[0], *objFunVal, objFunGrad[0]);
    } else {
        mock_c()->actualCall("testObj");

        *objFunVal = pow(parameters[0] + 1.0, 2) + 42.0;

//        printf("f: %f %f\n", parameters[0], *objFunVal);

    }
    return 0;

}

void QuadraticTestProblem::logOptimizerFinished(double optimalCost, const double *optimalParameters, double masterTime, int exitStatus)
{
    mock_c()->actualCall("logFinish")->withIntParameters("exitStatus", exitStatus);

    this->optimalCost = optimalCost;
    this->optimalParameter = optimalParameters[0];

//    printf("f(x) %f x %f t %f s %d\n", optimalCost, optimalParameters[0], masterTime, exitStatus);

}

QuadraticTestProblem::~QuadraticTestProblem()
{
    delete[] initialParameters;
    delete[] parametersMin;
    delete[] parametersMax;
}
