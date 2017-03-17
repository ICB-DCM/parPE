#include "optimizationProblem.h"
#include <stdlib.h>
#include <string.h>

OptimizationProblem *optimizationProblemNew()
{
    OptimizationProblem *problem = malloc(sizeof(*problem));
    memset(problem, 0, sizeof(*problem));

    return problem;
}

MultiStartOptimization *multiStartOptimizationNew()
{
    MultiStartOptimization *ms = malloc(sizeof(*ms));
    memset(ms, 0, sizeof(*ms));

    return ms;

}
