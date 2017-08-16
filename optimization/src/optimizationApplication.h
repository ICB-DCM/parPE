#ifndef OPTIMIZATIONAPPLICATION_H
#define OPTIMIZATIONAPPLICATION_H

#include "optimizationProblem.h"
#include "optimizationResultWriter.h"
#include <getopt.h>
#include <string>

/**
 * @brief The OptimizationApplication class parses command line arguments,
 * initializes MPI in required, opens data and results files and starts an
 * optimization
 */

class OptimizationApplication {
  public:
    OptimizationApplication();

    virtual int init(int argc, char **argv);

    virtual int parseOptions(int argc, char **argv);

    static void initMPI(int *argc, char ***argv);

    virtual void initProblem(const char *inFileArgument,
                             const char *outFileArgument) = 0;

    virtual void destroyProblem() {}

    virtual int run();

    virtual int runMaster() { return 0; }

    virtual void runWorker() {}

    virtual void runSingleMpiProcess() {}

    virtual void finalizeTiming(clock_t begin);

    ~OptimizationApplication();

    int getMpiRank();
    int getMpiCommSize();

    const char *dataFileName;
    char *resultFileName;
    OptimizationProblem *problem;
    OptimizationResultWriter *resultWriter;

    // command line option parsing
    const char *shortOptions = "dhvt:o:";
    struct option const longOptions[7] = {
        {"debug", no_argument, NULL, 'd'},
        {"print-worklist", no_argument, NULL, 'p'},
        {"help", no_argument, NULL, 'h'},
        {"version", no_argument, NULL, 'v'},
        {"task", required_argument, NULL, 't'},
        {"outfile-prefix", required_argument, NULL, 'o'},
        {NULL, 0, NULL, 0}};

    typedef enum operationType_tag {
        OP_TYPE_PARAMETER_ESTIMATION,
        OP_TYPE_GRADIENT_CHECK
    } operationTypeEnum;

    operationTypeEnum opType = OP_TYPE_PARAMETER_ESTIMATION;
};

#endif // OPTIMIZATIONAPPLICATION_H
