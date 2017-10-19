#include <stdio.h>

#include "optimizationOptions.h"
#include "optimizationProblem.h"
#include "steadystateProblem.h"

int main(int argc, char **argv) {
    ExampleSteadystateProblem problem = ExampleSteadystateProblem();

    printf("#########\n");
    printf("# IpOpt #\n");
    printf("#########\n");

    int status = parpe::getLocalOptimum(&problem);

    printf("#########\n");
    printf("# CERES #\n");
    printf("#########\n");

    problem.optimizationOptions->optimizer = parpe::OPTIMIZER_CERES;
    status += parpe::getLocalOptimum(&problem);

    return status;
}
