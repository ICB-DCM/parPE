#ifndef OPTIMIZATIONAPPLICATION_H
#define OPTIMIZATIONAPPLICATION_H

#include "optimizationProblem.h"
#include "optimizationResultWriter.h"
#include <string>

/**
 * @brief The OptimizationApplication class parses command line arguments, initializes MPI in required, opens data and results files and starts an optimization
 */

class OptimizationApplication
{
public:
    OptimizationApplication();
    OptimizationApplication(int argc, char **argv);

    static void initMPI(int *argc, char ***argv);

    virtual void initProblem(const char *inFileArgument, const char *outFileArgument) = 0;

    virtual void destroyProblem() {}

    virtual int run();

    virtual int runMaster() { return 0; }

    virtual void runWorker() {}

    ~OptimizationApplication();

    int getMpiRank();
    int getMpiCommSize();

    const char *dataFileName;
    const char *resultFileName;
    OptimizationProblem *problem;
    OptimizationResultWriter *resultWriter;


};

#endif // OPTIMIZATIONAPPLICATION_H
