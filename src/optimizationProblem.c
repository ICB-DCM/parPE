#include "optimizationProblem.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "misc.h"
#include "localOptimizationCeres.hpp"
#include "localOptimizationIpopt.h"

OptimizationProblem *optimizationProblemNew()
{
    OptimizationProblem *problem = malloc(sizeof(*problem));
    memset(problem, 0, sizeof(*problem));

    return problem;
}

/**
 * @brief getLocalOptimum
 * @param problem
 * @return int indicating status. 0: success, != 0: failure
 */

int getLocalOptimum(OptimizationProblem *problem)
{
    switch (problem->optimizer) {
    case OPTIMIZER_CERES:
        return getLocalOptimumCeres(problem);
    case OPTIMIZER_IPOPT:
        return getLocalOptimumIpopt(problem);
    default:
        abort();
    }
}

/**
 * @brief getLocalOptimumThreadWrapper wrapper for using getLocalOptimum with pThreads.
 * @param problem
 * @return Pointer to int indicating status. 0: success, != 0: failure
 */

void *getLocalOptimumThreadWrapper(void *optimizationProblemVp)
{
    OptimizationProblem *problem = (OptimizationProblem *) optimizationProblemVp;
    int *result = malloc(sizeof(*result));
    *result = getLocalOptimum(problem);
    return result;
}


void runOptimizationsParallel(const OptimizationProblem **problems, int numProblems) {
    runInParallelAndWaitForFinish(getLocalOptimumThreadWrapper, (void**)problems, numProblems);
}

void getRandomStartingpoint(const double *min, const double *max, int numParameters, double *buffer)
{
    fillArrayRandomDoubleIndividualInterval(min, max, numParameters, buffer);
}

void optimizationProblemGradientCheck(OptimizationProblem *problem, const int parameterIndices[], int numParameterIndices, double epsilon)
{
    double fc = 0; // f(theta)
    double *theta = problem->initialParameters;

    double *gradient = malloc(sizeof(double) * problem->numOptimizationParameters);
    problem->objectiveFunctionGradient(problem, theta, &fc, gradient);

    double *thetaTmp = malloc(sizeof(double) * problem->numOptimizationParameters);
    memcpy(thetaTmp, theta, sizeof(double) * problem->numOptimizationParameters);

    printf("Index\tGradient\tfd_f\t\t(delta)\t\tfd_c\t\t(delta)\t\tfd_b\t\t(delta)\n");

    for(int i = 0; i < numParameterIndices; ++i) {
        int curInd = parameterIndices[i];
        double fb = 0, ff = 0; // f(theta + eps) , f(theta - eps)

        thetaTmp[curInd] = theta[curInd] + epsilon;
        problem->objectiveFunction(problem, thetaTmp, &ff);

        thetaTmp[curInd] = theta[curInd] - epsilon;
        problem->objectiveFunction(problem, thetaTmp, &fb);

        double fd_f = (ff - fc) / epsilon;

        double fd_b = (fc - fb) / epsilon;

        double fd_c = (ff - fb) / (2 * epsilon);

        thetaTmp[curInd] = theta[curInd];

        printf("%d\tg: %f\tfd_f: %f\t(%f)\tfd_c: %f\t(%f)\tfd_b: %f\t(%f)",
               curInd, gradient[curInd],
               fd_f, gradient[curInd] - fd_f,
               fd_c, gradient[curInd] - fd_c,
               fd_b, gradient[curInd] - fd_b);
        printf("\t\tfb: %f\tfc: %f\tff: %f\t\n", fb, fc, ff);
    }

    free(gradient);
    free(thetaTmp);
}
