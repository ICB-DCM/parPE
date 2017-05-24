#ifndef QUADRATIC_TEST_PROBLEM_H
#define QUADRATIC_TEST_PROBLEM_H

#include "optimizationProblem.h"

class QuadraticTestProblem : public OptimizationProblem {
public:

    QuadraticTestProblem();
    /*
     * Test Problem for minimization
     * f (x) = (x+1)^2 + 42 = x^2  + 2x + 1 + 42
     * f'(x) = 2x + 2 = 0 <=> x = -1
     * f(-1) = 42
     */

    int evaluateObjectiveFunction(const double *parameters, double *objFunVal, double *objFunGrad);

    void logOptimizerFinished(double optimalCost,
                                const double *optimalParameters,
                                double masterTime,
                                int exitStatus);

    ~QuadraticTestProblem();


    double optimalCost;
    double optimalParameter;
};

#endif
