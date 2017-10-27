#include <cstdio>
#include "optimizationOptions.h"
#include "steadystateProblem.h"
#include <string>
#include <iostream>
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
    if(argc != 2) {
        std::cerr<<"Error: wrong number of arguments. Exactly one argument for data file expected.";
        return EXIT_FAILURE;
    }

    std::string dataFileName = argv[1];
    ExampleSteadystateProblem problem = ExampleSteadystateProblem(dataFileName);

    int status = 0;

    printf("#########\n");
    printf("# IpOpt #\n");
    printf("#########\n");

    status += parpe::getLocalOptimum(&problem);

    printf("#########\n");
    printf("# CERES #\n");
    printf("#########\n");

    parpe::OptimizationOptions options = problem.getOptimizationOptions();
    options.optimizer = parpe::OPTIMIZER_CERES;
    problem.setOptimizationOptions(options);

    status += parpe::getLocalOptimum(&problem);

#ifdef PARPE_DLIB_ENABLED
    printf("#########\n");
    printf("# Dlib  #\n");
    printf("#########\n");

    options.optimizer = parpe::OPTIMIZER_DLIB;
    problem.setOptimizationOptions(options);

    status += parpe::getLocalOptimum(&problem);
#endif

    return status;
}
