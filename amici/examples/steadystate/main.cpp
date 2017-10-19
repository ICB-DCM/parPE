#include <stdio.h>

#include "optimizationOptions.h"
#include "optimizationProblem.h"
#include "steadystateProblem.h"

int main(int argc, char **argv) {
    ExampleSteadystateProblem problem = ExampleSteadystateProblem();

    printf("#########\n");
    printf("# IpOpt #\n");
    printf("#########\n");

    int status = parPE::getLocalOptimum(&problem);

    printf("#########\n");
    printf("# CERES #\n");
    printf("#########\n");

    problem.optimizationOptions->optimizer = parPE::OPTIMIZER_CERES;
    status += parPE::getLocalOptimum(&problem);

    return status;
}
