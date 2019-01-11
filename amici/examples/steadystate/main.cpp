#include "steadystateProblem.h"
#include <optimizationOptions.h>
#include <string>
#include <iostream>
#include <cstdio>

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

    std::string dataFileName = argv[1]; // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
    ExampleSteadystateProblem problem(dataFileName);
    parpe::OptimizationOptions options = problem.getOptimizationOptions();

    int status = 0;

    printf("#########\n");
    printf("# IpOpt #\n");
    printf("#########\n");

    status += parpe::getLocalOptimum(&problem);

    printf("#########\n");
    printf("# CERES #\n");
    printf("#########\n");

    options.optimizer = parpe::optimizerName::OPTIMIZER_CERES;
    problem.setOptimizationOptions(options);

    status += parpe::getLocalOptimum(&problem);


#ifdef PARPE_DLIB_ENABLED
    printf("#########\n");
    printf("# Dlib  #\n");
    printf("#########\n");

    options.optimizer = parpe::optimizerName::OPTIMIZER_DLIB;
    problem.setOptimizationOptions(options);

    status += parpe::getLocalOptimum(&problem);
#endif

#ifdef PARPE_TOMS611_ENABLED
    printf("#########\n");
    printf("# TOMS611  #\n");
    printf("#########\n");

    options.optimizer = parpe::optimizerName::OPTIMIZER_TOMS611;
    problem.setOptimizationOptions(options);

    status += parpe::getLocalOptimum(&problem);
#endif

    return status;
}

//TODO


//void ExampleSteadystateProblem::logOptimizerFinished(
//    double optimalCost, const double *optimalParameters, double masterTime,
//    int exitStatus) {
//    printf("Minimal cost: %f\n", optimalCost);
//    printf("Optimal parameters  : ");
//    printArray(optimalParameters, model->np);
//    printf("\n");
//    printf("True parameters were: ");

//    hsize_t length;
//    double *ptrue;
//    AMI_HDF5_getDoubleArrayAttribute(fileId, "/parameters/", "ptrue", &ptrue,
//                                     &length);
//    printArray(ptrue, length);
//    delete ptrue;
//    printf("\n");

//    printf("Wall time (min): %f\n", masterTime / 60);
//}
