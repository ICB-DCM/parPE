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

    parpe::OptimizationOptions options = problem.getOptimizationOptions();
    options.optimizer = parpe::OPTIMIZER_CERES;
    problem.setOptimizationOptions(options);

    status += parpe::getLocalOptimum(&problem);

    return status;
}
