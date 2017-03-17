#include "optimizationProblem.h"
#include <stdlib.h>
#include <string.h>

#include "localOptimizationCeres.hpp"
#include "localOptimizationIpopt.h"

OptimizationProblem *optimizationProblemNew()
{
    OptimizationProblem *problem = malloc(sizeof(*problem));
    memset(problem, 0, sizeof(*problem));

    return problem;
}


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


void *getLocalOptimumThreadWrapper(void *optimizationProblemVp)
{
    OptimizationProblem *problem = (OptimizationProblem *) optimizationProblemVp;
    int *result = malloc(sizeof(*result));
    *result = getLocalOptimum(problem);
    return result;
}
