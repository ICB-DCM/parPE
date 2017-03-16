#include "CppUTest/TestHarness_c.h"
#include "CppUTestExt/MockSupport_c.h"

#include "localOptimizationIpopt.h"
#include "testingMisc.h"

double optimalCost;

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


TEST_GROUP_C_SETUP(localOptimizationIpopt) {
}

TEST_GROUP_C_TEARDOWN(localOptimizationIpopt) {
    mock_c()->checkExpectations();
    mock_c()->clear();
}

TEST_C(localOptimizationIpopt, testOptimization) {
    double theta0[] = {100};
    double thetaMin[] = {-1e5};
    double thetaMax[] = {1e5};

    OptimizationProblem *problem = optimizationProblemNew();

    problem->numOptimizationParameters = 1;
    problem->initialParameters = theta0;
    problem->parametersMin = thetaMin;
    problem->parametersMax = thetaMax;
    problem->objectiveFunction = testObj;
    problem->objectiveFunctionGradient = testObjGrad;
    problem->maxOptimizerIterations = 12;
    problem->logOptimizerFinished = logFinish;

    mock_c()->expectOneCall("logFinish")->withIntParameters("exitStatus", 0);
    getLocalOptimumIpopt(problem);

    // free(problem); //?
}


