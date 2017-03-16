#ifndef QUADRATIC_TEST_PROBLEM_H
#define QUADRATIC_TEST_PROBLEM_H

#include "optimizationProblem.h"

double optimalCost;

int testObj(void *problem, const double *parameters, double *result);

int testObjGrad(void *problem, const double *parameters, double *objFunVal, double *objFunGrad);

void logFinish(void *problem, double finalNegLogLikelihood, double masterTime, int exitStatus);

OptimizationProblem *getQuadraticTestProblem();


#endif
