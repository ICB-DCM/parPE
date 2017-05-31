#include "optimizationProblem.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "misc.h"
#include "localOptimizationCeres.h"
#include "localOptimizationIpopt.h"

/**
 * @brief getLocalOptimum
 * @param problem
 * @return int indicating status. 0: success, != 0: failure
 */

int getLocalOptimum(OptimizationProblem *problem)
{
    switch (problem->optimizationOptions->optimizer) {
    case OPTIMIZER_CERES:
        return getLocalOptimumCeres(problem);
    case OPTIMIZER_IPOPT:
        return getLocalOptimumIpopt(problem);
    default:
        return getLocalOptimumIpopt(problem);
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
    int *result = new int;
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

    double *gradient = new double[problem->numOptimizationParameters];
    problem->evaluateObjectiveFunction(theta, &fc, gradient);

    double *thetaTmp = new double[problem->numOptimizationParameters];
    memcpy(thetaTmp, theta, sizeof(double) * problem->numOptimizationParameters);

    printf("Index\tGradient\tfd_f\t\t(delta)\t\tfd_c\t\t(delta)\t\tfd_b\t\t(delta)\n");

    for(int i = 0; i < numParameterIndices; ++i) {
        int curInd = parameterIndices[i];
        double fb = 0, ff = 0; // f(theta + eps) , f(theta - eps)

        thetaTmp[curInd] = theta[curInd] + epsilon;
        problem->evaluateObjectiveFunction(thetaTmp, &ff, NULL);

        thetaTmp[curInd] = theta[curInd] - epsilon;
        problem->evaluateObjectiveFunction(thetaTmp, &fb, NULL);

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

    delete[] gradient;
    delete[] thetaTmp;
}

int OptimizationProblem::evaluateObjectiveFunction(const double *parameters, double *objFunVal, double *objFunGrad)
{
    return 0;
}

int OptimizationProblem::intermediateFunction(int alg_mod, int iter_count, double obj_value, double inf_pr, double inf_du, double mu, double d_norm, double regularization_size, double alpha_du, double alpha_pr, int ls_trials)
{
    return 0;
}

void OptimizationProblem::logObjectiveFunctionEvaluation(const double *parameters, double objectiveFunctionValue, const double *objectiveFunctionGradient, int numFunctionCalls, double timeElapsed)
{

}

void OptimizationProblem::logOptimizerFinished(double optimalCost, const double *optimalParameters, double masterTime, int exitStatus)
{

}

OptimizationProblem::OptimizationProblem()
{
    numOptimizationParameters = 0;
    initialParameters = NULL;
    parametersMin = NULL;
    parametersMax = NULL;

    optimizationOptions = NULL;
}

OptimizationOptions::OptimizationOptions()
{
    optimizer = OPTIMIZER_IPOPT;
    logFile = NULL;
    printToStdout = true;
    maxOptimizerIterations = 1000;
}

OptimizationOptions *OptimizationOptions::fromHDF5(const char *fileName)
{
    return NULL; // TODO
}
