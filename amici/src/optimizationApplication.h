#ifndef OPTIMIZATIONAPPLICATION_H
#define OPTIMIZATIONAPPLICATION_H

#include <loadBalancerMaster.h>
#include <loadBalancerWorker.h>

#include "multiConditionProblem.h"
#include "hierarchicalOptimization.h"

#include <getopt.h>
#include <string>

namespace parpe {

/**
 * @brief The OptimizationApplication class parses command line arguments,
 * initializes MPI if required, opens data and results files and starts an
 * optimization
 */
// TODO: DistributedOptimizationApplication
class OptimizationApplication {
  public:
    OptimizationApplication() = default;

    virtual ~OptimizationApplication();

    /**
     * @brief User-provided problem initialization.
     * Must set OptimizationApplication::problem, OptimizationApplication::multiStartOptimization and should set
     * OptimizationApplication::resultWriter
     * @param inFileArgument
     * @param outFileArgument
     */
    virtual void initProblem(std::string inFileArgument,
                             std::string outFileArgument) = 0;

    /**
     * @brief Start the optimization run. Must only be called once.
     * @return status code; 0 on success
     */
    int run(int argc, char **argv);

    /**
     * @brief This is run by the MPI rank 0 process when started with multiple
     * processes.
     * @return
     */
    virtual void runMaster();

    /**
     * @brief Code to be run on worker processes. Waits for jobs to be sent to
     * messageHandler()
     */
    virtual void runWorker();

    /**
     * @brief Code to be run if the application is running on only 1 process
     */
    virtual void runSingleMpiProcess();

  protected:

    /**
     * @brief Writes the total programm runtime
     * @param begin
     */
    virtual void finalizeTiming(double wallTimeSeconds, double cpuTimeSeconds);

    std::string
    processResultFilenameCommandLineArgument(const char *commandLineArg);

    bool isMaster();

    bool isWorker();

    /**
     * @brief Initialize MPI
     * @param argc
     * @param argv
     */
    static void initMPI(int *argc, char ***argv);

    /**
     * @brief Parse command line Options.
     * Must be called before any other functions.
     * Initializes MPI if not already done.
     * @param argc
     * @param argv
     * @return
     */
    virtual int parseOptions(int argc, char **argv);

    virtual void printUsage(char* const argZero);

    virtual void logParPEVersion(hid_t file_id) const;

private:
    /**
     * @brief initialize MPI, mutexes, ...
     * @param argc
     * @param argv
     */
    int init(int argc, char **argv);

    void runMultiStarts();

protected:
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

    // the need to be filled in by sub
    std::unique_ptr<MultiStartOptimizationProblem> multiStartOptimizationProblem;
    std::unique_ptr<OptimizationProblem> problem;
    hid_t file_id;
    operationTypeEnum opType = OP_TYPE_PARAMETER_ESTIMATION;
    LoadBalancerMaster loadBalancer;
};


/**
 * @brief CPU time for whole application run
 * @param timeInSeconds
 */
void saveTotalCpuTime(hid_t file_id, const double timeInSeconds);


} // namespace parpe

#endif // OPTIMIZATIONAPPLICATION_H
