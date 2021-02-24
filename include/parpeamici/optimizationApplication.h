#ifndef OPTIMIZATIONAPPLICATION_H
#define OPTIMIZATIONAPPLICATION_H

#include <parpecommon/parpeConfig.h>
#include <parpeamici/multiConditionProblem.h>
#include <parpeamici/hierarchicalOptimization.h>

#ifdef PARPE_ENABLE_MPI
#include <parpeloadbalancer/loadBalancerMaster.h>
#include <parpeloadbalancer/loadBalancerWorker.h>
#endif

#include <getopt.h>
#include <string>

namespace parpe {

#ifndef PARPE_ENABLE_MPI
// Workaround to allow building without MPI. Should be cleaned up.
using LoadBalancerMaster = int;
using LoadBalancerWorker = int;
#endif

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
     * Must set OptimizationApplication::problem,
     * OptimizationApplication::multiStartOptimization and should set
     * OptimizationApplication::resultWriter
     * @param inFileArgument
     * @param outFileArgument
     */
    virtual void initProblem(std::string inFileArgument,
                             std::string outFileArgument) = 0;

    /**
     * @brief Start the optimization run. Must only be called once.
     * Initializes MPI if not already done.
     * Must be called before any other functions.
     * @return status code; 0 on success
     */
    int run(int argc, char **argv);

    /**
     * @brief This is run by the MPI rank 0 process when started with multiple
     * processes.
     */
    virtual void runMaster();

#ifdef PARPE_ENABLE_MPI
    /**
     * @brief Code to be run on worker processes. Waits for jobs to be sent to
     * messageHandler()
     */
    virtual void runWorker();
#endif

    /**
     * @brief Code to be run if the application is running on only 1 process
     */
    virtual void runSingleProcess();

  protected:

    /**
     * @brief Receives and writes the total programm runtime
     * @param begin
     */
    virtual void finalizeTiming(double wallTimeSeconds, double cpuTimeSeconds);

    /**
     * @brief processResultFilenameCommandLineArgument
     * @param commandLineArg
     * @return Result file name
     */
    std::string
    processResultFilenameCommandLineArgument(const char *commandLineArg);

    /**
     * @brief Are we running with MPI and are we master process?
     * @return
     */
    bool isMaster();

    /**
     * @brief Are we running with MPI and are we a worker?
     * @return
     */
    bool isWorker();

    /**
     * @brief Initialize MPI
     * @param argc
     * @param argv
     */
    static void initMPI(int *argc, char ***argv);

    /**
     * @brief Parse command line options before MPI_INIT is potentially called.
     *
     * Used e.g. to print usage information without first initialization MPI.
     *
     * Argv may contain extra MPI arguments.
     * @param argc
     * @param argv
     * @return
     */
    virtual int parseCliOptionsPreMpiInit(int argc, char **argv);

    /**
     * @brief Parse command line options after MPI_Init is called.
     *
     * Any MPI-related CLI arguments will be removed here.
     * @param argc
     * @param argv
     * @return
     */
    virtual int parseCliOptionsPostMpiInit(int argc, char **argv);

    /**
     * @brief Print CLI usage
     * @param argZero
     */
    virtual void printUsage(char* const argZero);

    virtual void logParPEVersion(const H5::H5File &file) const;

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
    const char *shortOptions = "dhvmt:o:s:";
    struct option const longOptions[10] = {
        {"debug", no_argument, NULL, 'd'},
        {"print-worklist", no_argument, NULL, 'p'},
        {"help", no_argument, NULL, 'h'},
        {"version", no_argument, NULL, 'v'},
        {"mpi", no_argument, NULL, 'm'},
        {"task", required_argument, NULL, 't'},
        {"gradient-check", required_argument, NULL, 'g'},
        {"outfile-prefix", required_argument, NULL, 'o'},
        {"first-start-idx", required_argument, NULL, 's'},
        {NULL, 0, NULL, 0}};

    enum class OperationType {
        parameterEstimation,
        gradientCheck
    };

    std::string dataFileName;
    std::string resultFileName;

    int first_start_idx {0};

    // the need to be filled in by sub
    std::unique_ptr<MultiStartOptimizationProblem> multiStartOptimizationProblem;
    std::unique_ptr<OptimizationProblem> problem;
    H5::H5File h5File = 0;
    OperationType operationType = OperationType::parameterEstimation;
    int num_parameter_checks = 1;
    LoadBalancerMaster loadBalancer;
    bool withMPI = false;
};


/**
 * @brief CPU time for whole application run
 * @param file
 * @param timeInSeconds
 */
void saveTotalCpuTime(const H5::H5File &file, const double timeInSeconds);


} // namespace parpe

#endif // OPTIMIZATIONAPPLICATION_H
