#include "optimizationApplication.h"
#include "LoadBalancerMaster.h"
#include "hdf5Misc.h"
#include "logging.h"
#include "misc.h"
#include <cstring>
#include <ctime>
#include <mpi.h>
#include <optimizationOptions.h>
#include <pthread.h>
#include <amici.h>
#include <amiciMisc.h>
#include <numeric>
#include <algorithm>
#include <random>

namespace parpe {

OptimizationApplication::OptimizationApplication()
    : OptimizationApplication(0, nullptr) {}

OptimizationApplication::OptimizationApplication(int argc, char **argv) {
    // TODO: check if initialized already
    initMPI(&argc, &argv);

    printMPIInfo();

    initHDF5Mutex();

    parseOptions(argc, argv);

    // Seed random number generator
    //    srand(1337);
    unsigned int seed = time(NULL);
    logmessage(LOGLVL_DEBUG, "Seeding RNG with %u", seed);
    srand(seed); // TODO to CLI

    amici::errMsgIdAndTxt = printAmiciErrMsgIdAndTxt;
    amici::warnMsgIdAndTxt = printAmiciWarnMsgIdAndTxt;
}

int OptimizationApplication::parseOptions(int argc, char **argv) {
    int c;

    while (1) {
        int optionIndex = 0;
        c = getopt_long(argc, argv, shortOptions, longOptions, &optionIndex);

        if (c == -1)
            break;

        switch (c) {
        case 'd':
            printDebugInfoAndWait();
            break;
        case 't':
            if (strcmp(optarg, "gradient_check") == 0)
                opType = OP_TYPE_GRADIENT_CHECK;
            break;
        case 'o':
            resultFileName = processResultFilenameCommandLineArgument(optarg);
            break;
        case 'v':
            printf("Version: %s\n", GIT_VERSION);
            return 1;
        case 'h':
            printUsage(argv[0]);
            return 1;
        default:
            printf("Unrecognized option: %c\n", c);
        }
    }

    if (optind < argc) {
        dataFileName = argv[optind++];
    } else {
        logmessage(LOGLVL_CRITICAL,
                   "Must provide input file as first and only argument to %s.",
                   argv[0]);
        return 1;
    }

    return 0;
}

void OptimizationApplication::printUsage(char * const argZero)
{
    printf("Usage: %s [OPTION]... FILE\n", argZero);
    printf("FILE: HDF5 data file");
    printf("Options: \n"
           "  -o, --outfile-prefix Prefix for result files (path + "
           "filename)\n"
           "  -t, --task    What to do? Parameter estimation (default) "
           "or check gradient ('gradient_check')\n"
           "  -h, --help    Print this help text\n"
           "  -v, --version Print version info\n");
}

void OptimizationApplication::initMPI(int *argc, char ***argv) {
    int mpiErr = MPI_Init(argc, argv);
    if (mpiErr != MPI_SUCCESS) {
        logmessage(LOGLVL_CRITICAL, "Problem initializing MPI. Exiting.");
        exit(1);
    }

    int mpiRank = getMpiRank();

    if (mpiRank == 0) {
        int commSize = getMpiCommSize();
        logmessage(LOGLVL_INFO, "Running with %d MPI processes.", commSize);
    }
}

int OptimizationApplication::run() {
    clock_t begin = clock();

    int status = 0;

    if (!dataFileName.size()) {
        logmessage(LOGLVL_CRITICAL,
                   "No input file provided. Must provide input file as first "
                   "and only argument or set "
                   "OptimizationApplication::inputFileName manually.");
        return 1;
    }
    initProblem(dataFileName, resultFileName);

    int commSize = getMpiCommSize();

    if (commSize > 1) {
        if (getMpiRank() == 0) {
            loadBalancer.run();
            status = runMaster();

            loadBalancer.terminate();
            loadBalancer.sendTerminationSignalToAllWorkers();
            finalizeTiming(begin);
            logmessage(LOGLVL_INFO, "Sent termination signal to workers.");

        } else {
            status = runWorker();
            finalizeTiming(begin);
        }
    } else {
        runSingleMpiProcess();

        finalizeTiming(begin);
    }

    return status;
}

int OptimizationApplication::runMaster() {
    switch (opType) {
    case OP_TYPE_GRADIENT_CHECK: {
        // startObjectiveFunctionGradientCheck(&problem);

        const int numParameterIndicesToCheck = 50;
        const double epsilon = 1e-6;

        // choose random parameters to check
        std::vector<int> parameterIndices(problem->getNumOptimizationParameters());
        std::iota(parameterIndices.begin(), parameterIndices.end(), 0);
        std::random_device rd;
        std::mt19937 g(rd());
        std::shuffle(parameterIndices.begin(), parameterIndices.end(), g);

        optimizationProblemGradientCheck(problem.get(), parameterIndices.data(),
                                         numParameterIndicesToCheck, epsilon);
        break;
    }
    case OP_TYPE_PARAMETER_ESTIMATION:
    default:
        // startParameterEstimation(&dataProvider);

        // if numStarts > 1: need to use multiple MPI
        // workers, otherwise simulation crashes due
        // to CVODES threading issues

        if (problem->getOptimizationOptions().numStarts > 0) {
            MultiConditionProblemMultiStartOptimization ms(
                problem->getOptimizationOptions().numStarts,
                problem->getOptimizationOptions().retryOptimization);
            ms.options = problem->getOptimizationOptions();
            ms.resultWriter = problem->resultWriter.get();
            ms.dp = problem->getDataProvider();
            ms.loadBalancer = &loadBalancer;
            ms.run();
        }

        break;
    }

    return 0;
}

int OptimizationApplication::runSingleMpiProcess() {
    // TODO: also for gradientCheck
    // run serially
    if (problem->getOptimizationOptions().numStarts > 0) {
        MultiConditionProblemMultiStartOptimization ms(
            problem->getOptimizationOptions().numStarts,
            problem->getOptimizationOptions().retryOptimization);
        ms.options = problem->getOptimizationOptions();
        ms.resultWriter = problem->resultWriter.get();
        ms.dp = problem->getDataProvider();
        return ms.run();
    } else {
        return getLocalOptimum(problem.get());
    }
}

void OptimizationApplication::finalizeTiming(clock_t begin) {
    // wall-time for current process
    clock_t end = clock();
    double wallTimeSeconds = (double)(end - begin) / CLOCKS_PER_SEC;

    // total run-time
    double totalTimeInSeconds = 0;
    MPI_Reduce(&wallTimeSeconds, &totalTimeInSeconds, 1, MPI_DOUBLE, MPI_SUM, 0,
               MPI_COMM_WORLD);

    int mpiRank = getMpiRank();

    if (mpiRank == 0) {
        logmessage(LOGLVL_INFO, "Walltime: %fs, total compute time:%fs",
                   wallTimeSeconds, totalTimeInSeconds);
        if (resultWriter)
            resultWriter->saveTotalCpuTime(totalTimeInSeconds);
    }
}

std::string OptimizationApplication::processResultFilenameCommandLineArgument(
    const char *commandLineArg) {
    std::size_t bufSize = 1024;
    char tmpFileName[bufSize];
    snprintf(tmpFileName, bufSize, "%s_rank%05d.h5", commandLineArg,
             getMpiRank());
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
    problem.reset(nullptr);
    resultWriter.reset(nullptr);
    MPI_Finalize();
}

} // namespace parpe
