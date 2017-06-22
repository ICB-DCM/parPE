#ifndef OPTIMIZATIONAPPLICATION_H
#define OPTIMIZATIONAPPLICATION_H

#include "optimizationProblem.h"
#include <string>

/**
 * @brief The OptimizationApplication class parses command line arguments, initializes MPI in required, opens data and results files and starts an optimization
 */

class OptimizationApplication
{
public:
    OptimizationApplication();
    OptimizationApplication(OptimizationProblem *problem, int argc, char **argv);

    static void initMPI(int *argc, char ***argv);

    int run();

    ~OptimizationApplication();

    const char *dataFileName;
    const char *resultFileName;
    OptimizationProblem *problem;


};

#endif // OPTIMIZATIONAPPLICATION_H
