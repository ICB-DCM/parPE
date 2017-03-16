#include "CppUTest/TestHarness_c.h"
#include "CppUTestExt/MockSupport_c.h"

#include "quadraticTestProblem.h"

double theta0[] = {100};
double thetaMin[] = {-1e5};
double thetaMax[] = {1e5};

int testObj(void *problem, const double *parameters, double *result){
    *result = parameters[0] * parameters[0] + 42;
    return 0;
}

int testObjGrad(void *problem, const double *parameters, double *objFunVal, double *objFunGrad) {
    objFunGrad[0] = 2 * parameters[0];
    return 0;
}

void logFinish(void *problem, double finalNegLogLikelihood, double masterTime, int exitStatus) {
    optimalCost = finalNegLogLikelihood;
    // printf("%f %f %d\n", finalNegLogLikelihood, masterTime, exitStatus);
    mock_c()->actualCall("logFinish")->withIntParameters("exitStatus", exitStatus);
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
