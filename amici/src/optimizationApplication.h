#ifndef OPTIMIZATIONAPPLICATION_H
#define OPTIMIZATIONAPPLICATION_H

#include "multiConditionProblem.h"
#include "multiConditionProblemResultWriter.h"
#include <LoadBalancerMaster.h>
#include <LoadBalancerWorker.h>
#include <getopt.h>
#include <string>
/**
 * @brief The OptimizationApplication class parses command line arguments,
 * initializes MPI in required, opens data and results files and starts an
 * optimization
 */
// TODO: DistributedOptimizationApplication
class OptimizationApplication {
  public:
    OptimizationApplication();

    OptimizationApplication(int argc, char **argv);

    ~OptimizationApplication();

    /**
     * @brief User-provided problem initialization
     * @param inFileArgument
     * @param outFileArgument
     */
    virtual void initProblem(std::string inFileArgument,
                             std::string outFileArgument) = 0;

    /**
     * @brief Start the optimization run
     * @return
     */
    int run();

    /**
     * @brief This is run by the MPI rank 0 process when started with multiple
     * processes.
     * @return
     */
    virtual int runMaster();

    /**
     * @brief Code to be run on worker processes. Waits for jobs to be sent to
     * messageHandler()
     */
    virtual int runWorker() {
        // TODO: Move out of here
        problem->run();

        return 0;
    }

    /**
     * @brief Code to be run if the application is running on only 1 process
     */
    virtual int runSingleMpiProcess();

    /**
     * @brief Writes the total programm runtime
     * @param begin
     */
    virtual void finalizeTiming(clock_t begin);

    std::string
    processResultFilenameCommandLineArgument(const char *commandLineArg);

    bool isMaster();

    bool isWorker();

    static int getMpiRank();

    static int getMpiCommSize();

  protected:
    /**
     * @brief Initialize MPI
     * @param argc
     * @param argv
     */
    static void initMPI(int *argc, char ***argv);

    /**
     * @brief Parse command line Options.
     * Must be called before any other functions. Initializes MPI if not already
     * done.
     * @param argc
     * @param argv
     * @return
     */
    virtual int parseOptions(int argc, char **argv);

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

    std::string dataFileName;
    std::string resultFileName;
    MultiConditionProblem *problem = nullptr;
    MultiConditionProblemResultWriter *resultWriter = nullptr;
    operationTypeEnum opType = OP_TYPE_PARAMETER_ESTIMATION;
    LoadBalancerMaster loadBalancer;
};

#endif // OPTIMIZATIONAPPLICATION_H
