#ifndef QUADRATIC_TEST_PROBLEM_H
#define QUADRATIC_TEST_PROBLEM_H

#include "optimizationProblem.h"

extern double optimalCost, optimalParameter;

int testObj(void *problem, const double *parameters, double *result);

int testObjGrad(void *problem, const double *parameters, double *objFunVal, double *objFunGrad);

void logFinish(void *problem, double optimalCost, const double *optimalParameters, double masterTime, int exitStatus);

OptimizationProblem *getQuadraticTestProblem();


#endif
