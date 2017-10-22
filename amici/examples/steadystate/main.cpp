#include <stdio.h>
#include "optimizationOptions.h"
#include "steadystateProblem.h"

/**
 * @file
 *
 * This is an example for parameter estimation for the Steadystate ODE example model included in AMICI.
 * It demonstrates how to use IpOpt or CERES to solve a ODE-constrained optimization problem for which
 * the ODE system has been implemented in AMICI.
 * For cases where the ODE has to be evaluated several times per objective function evaluation
 * see examples example_steadystate_parallel and example_steadystate_multicondition.
 */

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
