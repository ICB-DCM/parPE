#include <parpeamici/optimizationApplication.h>

#include <parpecommon/hdf5Misc.h>
#include <parpecommon/logging.h>
#include <parpecommon/misc.h>
#include <parpeoptimization/optimizationOptions.h>
#include <parpecommon/parpeVersion.h>
#include <parpeamici/amiciMisc.h>

#ifdef PARPE_ENABLE_MPI
#include <mpi.h>
#endif

#include <cstring>
#include <ctime>
#include <numeric>
#include <algorithm>
#include <random>
#include <csignal>
#include <cstdlib>

namespace parpe {

struct sigaction act;
struct sigaction oldact;

void signalHandler(int sig) {
    logmessage(loglevel::critical, "Caught signal %d ", sig);
    printBacktrace();

    // restore previous
    oldact.sa_flags = SA_SIGINFO;
    (*oldact.sa_sigaction)(sig, nullptr, nullptr);
}

int OptimizationApplication::init(int argc, char **argv) {
    // reduce verbosity
    if(std::getenv("PARPE_NO_DEBUG"))
        minimumLogLevel = loglevel::info;

    if(auto status = parseCliOptionsPreMpiInit(argc, argv))
        return status;

    // install signal handler for backtrace on error
    sigaction(SIGSEGV, &act, &oldact);
    sigaction(SIGHUP, &act, nullptr);

    if(withMPI && !getMpiActive())
        initMPI(&argc, &argv);

    printMPIInfo();
    initHDF5Mutex();

    return parseCliOptionsPostMpiInit(argc, argv);
}

void OptimizationApplication::runMultiStarts() const
{
    // TODO: use uniqe_ptr, not ref
    MultiStartOptimization optimizer(*multiStartOptimizationProblem, true,
                                     first_start_idx);
    optimizer.run();
}

int OptimizationApplication::parseCliOptionsPreMpiInit(int argc, char **argv)
{
    while (true) {
        int optionIndex = 0;
        auto c = getopt_long(argc, argv, shortOptions, longOptions, &optionIndex);

        if (c == -1)
            break; // no more options

        switch (c) {
        case 'm':
            withMPI = true;
            break;
        case 'h':
            printUsage(argv[0]);
            return 1;
        }
    }
    return 0;
}

int OptimizationApplication::parseCliOptionsPostMpiInit(int argc, char **argv) {
    // restart from first argument
    optind = 1;

    while (true) {
        int optionIndex = 0;
        int c = getopt_long(argc, argv, shortOptions, longOptions, &optionIndex);

        if (c == -1)
            break; // no more options

        switch (c) {
        case 'd':
            printDebugInfoAndWait();
            break;
        case 't':
            if (strcmp(optarg, "gradient_check") == 0)
                operationType = OperationType::gradientCheck;
            break;
        case 'o':
            resultFileName = processResultFilenameCommandLineArgument(optarg);
            break;
        case 's':
            first_start_idx = atoi(optarg);
            break;
        case 'v':
            printf("Version: %s\n", PARPE_VERSION);
            return 1;
        case 'h':
        case 'm':
            continue;
        default:
            printf("Unrecognized option: %c\n", c);
            exit(EXIT_FAILURE);
        }
    }

    if (optind < argc) {
        dataFileName = argv[optind++];
    } else {
        logmessage(loglevel::critical,
                   "Must provide input file as first and only argument to %s.",
                   argv[0]);
        return 1;
    }

    return 0;
}

void OptimizationApplication::printUsage(char * const argZero)
{
    printf("Usage: %s [OPTION]... FILE\n\n", argZero);
    printf("FILE: HDF5 data file\n\n");
    printf("Options: \n"
           "  -o, --outfile-prefix  Prefix for result files (path + "
           "filename)\n"
           "  -t, --task            What to do? Parameter estimation (default) "
           "or check gradient ('gradient_check')\n"
           "  -s, --first-start-idx Starting point index for first optimization\n"
           "  -m, --mpi             Enable MPI (default: off)\n"
           "  -h, --help            Print this help text\n"
           "  -v, --version         Print version info\n");
    printf("\nSupported optimizers:\n");
    printAvailableOptimizers("  ");
}

void OptimizationApplication::logParPEVersion(H5::H5File const& file) const
{
    hdf5WriteStringAttribute(file, "/", "PARPE_VERSION", PARPE_VERSION);
}

void OptimizationApplication::initMPI(int *argc, char ***argv) {
#ifdef PARPE_ENABLE_MPI
    int thread_support_provided = 0;
    int mpiErr = MPI_Init_thread(argc, argv, MPI_THREAD_MULTIPLE,
                                 &thread_support_provided);

    if (thread_support_provided != MPI_THREAD_MULTIPLE)
        throw std::runtime_error("MPI_THREAD_MULTIPLE not supported?");

    if (mpiErr != MPI_SUCCESS) {
        logmessage(loglevel::critical, "Problem initializing MPI. Exiting.");
        exit(1);
    }

    int mpiRank = getMpiRank();

    if (mpiRank == 0) {
        int commSize = getMpiCommSize();
        logmessage(loglevel::info, "Running with %d MPI processes.", commSize);
    }
#endif
}

int OptimizationApplication::run(int argc, char **argv) {
    // start Timers
    WallTimer wallTimer;
    CpuTimer cpuTimer;

    int status = init(argc, argv);
    if(status)
        return status;

    if (dataFileName.empty()) {
        logmessage(loglevel::critical,
                   "No input file provided. Must provide input file as first "
                   "and only argument or set "
                   "OptimizationApplication::inputFileName manually.");
        return 1;
    }

    initProblem(dataFileName, resultFileName);

#ifdef PARPE_ENABLE_MPI
    if (getMpiCommSize() > 1) {
        if (getMpiRank() == 0) {
            loadBalancer.run();
            runMaster();

            loadBalancer.terminate();
            loadBalancer.sendTerminationSignalToAllWorkers();
            finalizeTiming(wallTimer.getTotal(), cpuTimer.getTotal());
            logmessage(loglevel::info, "Sent termination signal to workers.");

        } else {
            runWorker();
            finalizeTiming(wallTimer.getTotal(), cpuTimer.getTotal());
        }
    } else {
#endif
        runSingleProcess();

        finalizeTiming(wallTimer.getTotal(), cpuTimer.getTotal());
#ifdef PARPE_ENABLE_MPI
    }
#endif
    return status;
}

void OptimizationApplication::runMaster() {
    switch (operationType) {
    case OperationType::gradientCheck: {
        const int numParameterIndicesToCheck = 10000;
        optimizationProblemGradientCheckMultiEps(
                    problem.get(), numParameterIndicesToCheck);
        break;
    }
    case OperationType::parameterEstimation:
    default:
        runMultiStarts();
    }
}

#ifdef PARPE_ENABLE_MPI
void OptimizationApplication::runWorker() {
    // TODO: Move out of here
    LoadBalancerWorker lbw;
    lbw.run([this](std::vector<char> &buffer, int jobId) {
        // TODO: this is so damn ugly
        auto sgf = dynamic_cast<SummedGradientFunctionGradientFunctionAdapter<int>*>(problem->cost_fun_.get());
        if(sgf) {
            // non-hierarchical
            dynamic_cast<AmiciSummedGradientFunction*>(sgf->getWrappedFunction())->messageHandler(buffer, jobId);
        } else {
            // hierarchical
            auto hierarch = dynamic_cast<HierarchicalOptimizationWrapper *>(problem->cost_fun_.get());
            Expects(hierarch != nullptr);
            hierarch->getWrappedFunction()->messageHandler(buffer, jobId);
        }
    });
}
#endif

void OptimizationApplication::runSingleProcess() {
    // run serially
    switch (operationType) {
    case OperationType::gradientCheck: {
        const int numParameterIndicesToCheck = 10000;
        optimizationProblemGradientCheckMultiEps(
                    problem.get(), numParameterIndicesToCheck);
        break;
    }
    case OperationType::parameterEstimation:
    default:
        if (problem->getOptimizationOptions().numStarts > 0) {
            runMultiStarts();
        } else {
            getLocalOptimum(problem.get());
        }
    }
}

void OptimizationApplication::finalizeTiming(double wallTimeSeconds, double cpuTimeSeconds) {
#ifdef PARPE_ENABLE_MPI
    // wall-time for current MPI process
    // total run-time
    double totalCpuTimeInSeconds = 0;

    if(getMpiActive()) {
        MPI_Reduce(&cpuTimeSeconds, &totalCpuTimeInSeconds, 1, MPI_DOUBLE, MPI_SUM, 0,
                   MPI_COMM_WORLD);
    } else {
        totalCpuTimeInSeconds = cpuTimeSeconds;
    }

    if (getMpiRank() < 1) {
        logmessage(loglevel::info, "Walltime on master: %fs, CPU time of all processes: %fs",
                   wallTimeSeconds, totalCpuTimeInSeconds);
        saveTotalCpuTime(h5File, totalCpuTimeInSeconds);
    }
#else
    logmessage(LOGLVL_INFO, "Total walltime: %fs, CPU time: %fs",
               wallTimeSeconds, cpuTimeSeconds);
    saveTotalCpuTime(file_id, cpuTimeSeconds);
#endif
}

std::string OptimizationApplication::processResultFilenameCommandLineArgument(
        const char *commandLineArg) {
    std::size_t bufSize = 1024;
    char tmpFileName[bufSize];
    int rank = std::max(getMpiRank(), 0);
    snprintf(tmpFileName, bufSize, "%s_rank%05d.h5", commandLineArg, rank);
    return tmpFileName;

    /*
        // create directory for each compute node
        char procName[MPI_MAX_PROCESSOR_NAME];
        int procNameLen;
        MPI_Get_processor_name(procName, &procNameLen);

        createDirectoryIfNotExists(procName);

        char tmpFileName[200];
        strFormatCurrentLocaltime(tmpFileName, 200,
       "%%s/%Y-%m-%d_%H%M%S_%%05d.h5");

        int mpiRank;
        MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

        char *fileName = new char[MPI_MAX_PROCESSOR_NAME + 32];
        sprintf(fileName, tmpFileName, procName, mpiRank);
        std::string strFileName(fileName);
        delete[] fileName;

        return strFileName;
    */
}

bool OptimizationApplication::isMaster() { return getMpiRank() == 0; }

bool OptimizationApplication::isWorker() { return getMpiRank() > 0; }

OptimizationApplication::~OptimizationApplication() {
    // objects must be destroyed before MPI_Finalize is called
    // and Hdf5 mutex is destroyed
    [[maybe_unused]] auto lock = hdf5MutexGetLock();
    h5File.close();
    problem.reset(nullptr);
#ifdef PARPE_ENABLE_MPI
    if(withMPI && getMpiActive())
        MPI_Finalize();
#endif
}

void saveTotalCpuTime(H5::H5File const& file, const double timeInSeconds)
{
    hsize_t dims[1] = {1};

    [[maybe_unused]] auto lock = hdf5MutexGetLock();

    //std::string pathStr = rootPath + "/totalTimeInSec";
    std::string pathStr = "/totalTimeInSec";
    H5LTmake_dataset(file.getId(), pathStr.c_str(), 1, dims, H5T_NATIVE_DOUBLE,
                     &timeInSeconds);

}

} // namespace parpe
