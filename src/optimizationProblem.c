#include "optimizationProblem.h"
#include <stdlib.h>
#include <string.h>

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
