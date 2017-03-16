#include "CppUTest/TestHarness_c.h"
#include "CppUTestExt/MockSupport_c.h"

#include "quadraticTestProblem.h"
#include <math.h>
#include <stdio.h>

double quadraticTestProblemOptimalCost, quadraticTestProblemOptimalParameter;

double theta0[] = {-100};
double thetaMin[] = {-1e5};
double thetaMax[] = {1e5};

/*
 * Test Problem for minimization
 * f (x) = (x+1)^2 + 42 = x^2  + 2x + 1 + 42
 * f'(x) = 2x + 2 = 0 <=> x = -1
 * f(-1) = 42
 */

int testObj(void *problem, const double *parameters, double *objFunVal){

    mock_c()->actualCall("testObj");

    *objFunVal = pow(parameters[0] + 1.0, 2) + 42.0;

//    printf("f: %f %f\n", parameters[0], *objFunVal);

    return 0;
}

int testObjGrad(void *problem, const double *parameters, double *objFunVal, double *objFunGrad) {

    mock_c()->actualCall("testObjGrad");

    objFunVal[0] = pow(parameters[0] + 1.0, 2) + 42.0;
    objFunGrad[0] = 2.0 * parameters[0] + 2.0;

//    printf("g: %f %f %f\n", parameters[0], *objFunVal, objFunGrad[0]);

    return 0;
}

void logFinish(void *problem, double optimalCost, const double *optimalParameters, double masterTime, int exitStatus) {
    mock_c()->actualCall("logFinish")->withIntParameters("exitStatus", exitStatus);

    quadraticTestProblemOptimalCost = optimalCost;
    quadraticTestProblemOptimalParameter = optimalParameters[0];

//    printf("f(x) %f x %f t %f s %d\n", optimalCost, optimalParameters[0], masterTime, exitStatus);
}

OptimizationProblem *getQuadraticTestProblem()
{
    OptimizationProblem *problem = optimizationProblemNew();

    problem->numOptimizationParameters = 1;
    problem->initialParameters = theta0;
    problem->parametersMin = thetaMin;
    problem->parametersMax = thetaMax;
    problem->objectiveFunction = testObj;
    problem->objectiveFunctionGradient = testObjGrad;
    problem->maxOptimizerIterations = 12;
    problem->logOptimizerFinished = logFinish;

    return problem;
}
