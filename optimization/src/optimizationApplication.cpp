#include "optimizationApplication.h"
#include "logging.h"
#include "hdf5Misc.h"
#include "loadBalancerMaster.h"
#include <cstring>
#include <pthread.h>
#include <mpi.h>
#include <ctime>

OptimizationApplication::OptimizationApplication() : dataFileName(NULL), resultFileName(NULL), problem(NULL), resultWriter(NULL)
{

}

int OptimizationApplication::init(int argc, char **argv)
{
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

    return 0;
}

int OptimizationApplication::parseOptions(int argc, char **argv)
{
    int c;

    while (1) {
        int optionIndex = 0;
        c = getopt_long (argc, argv, shortOptions, longOptions, &optionIndex);

        if (c == -1)
            break;

        switch (c) {
        case 'd':
            printDebugInfoAndWait();
            break;
        case 't':
            if(strcmp(optarg, "gradient_check") == 0)
                opType = OP_TYPE_GRADIENT_CHECK;
            break;
        case 'o':
            resultFileName = new char[strlen(optarg) + 20];
            sprintf(resultFileName, "%s_rank%05d.h5", optarg, getMpiRank());
            break;
        case 'v':
            printf("Version: %s\n", GIT_VERSION);
            return 1;
        case 'h':
            printf("Usage: %s [OPTION]... FILE\n", argv[0]);
            printf("FILE: HDF5 data file");
            printf("Options: \n"
                   "  -o, --outfile-prefix Prefix for result files (path + filename)\n"
                   "  -t, --task    What to do? Parameter estimation (default) or check gradient ('gradient_check')\n"
                   "  -h, --help    Print this help text\n"
                   "  -v, --version Print version info\n"
                   );
            return 1;
        default:
            printf("Unrecognized option: %c\n", c);
        }
    }

    if(optind < argc) {
        dataFileName = argv[optind++];
    } else {
        logmessage(LOGLVL_CRITICAL, "Must provide input file as first and only argument to %s.", argv[0]);
        return 1;
    }

    return 0;
}

void OptimizationApplication::initMPI(int *argc, char ***argv)
{
    int mpiErr = MPI_Init(argc, argv);
    if(mpiErr != MPI_SUCCESS) {
        logmessage(LOGLVL_CRITICAL, "Problem initializing MPI. Exiting.");
        exit(1);
    }

    int mpiRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    if(mpiRank == 0) {
        int commSize;
        MPI_Comm_size(MPI_COMM_WORLD, &commSize);

        logmessage(LOGLVL_INFO, "Running with %d MPI processes.", commSize);
    }
}


int OptimizationApplication::run()
{
    clock_t begin = clock();

    int status = 0;

    if(!dataFileName) {
        logmessage(LOGLVL_CRITICAL, "No input file provided. Must provide input file as first and only argument or set OptimizationApplication::inputFileName manually.");
        return 1;
    }
    initProblem(dataFileName, resultFileName);

    int commSize;
    MPI_Comm_size(MPI_COMM_WORLD, &commSize);

    if(commSize > 1) {
        if(getMpiRank() == 0) {
            loadBalancerStartMaster();

            status = runMaster();

            loadBalancerTerminate();
            sendTerminationSignalToAllWorkers();
            finalizeTiming(begin);
            logmessage(LOGLVL_INFO, "Sent termination signal to workers.");

        } else {
            runWorker();
            finalizeTiming(begin);
        }
    } else {
        runSingleMpiProcess();

        finalizeTiming(begin);
    }

    return status;
}

void OptimizationApplication::finalizeTiming(clock_t begin)
{
    // wall-time for current process
    clock_t end = clock();
    double wallTimeSeconds = (double)(end - begin) / CLOCKS_PER_SEC;

    // total run-time
    double totalTimeInSeconds = 0;
    MPI_Reduce(&wallTimeSeconds, &totalTimeInSeconds, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    int mpiRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

    if(mpiRank == 0) {
        logmessage(LOGLVL_INFO, "Walltime: %fs, total compute time:%fs", wallTimeSeconds, totalTimeInSeconds);
        if(resultWriter)
            resultWriter->saveTotalCpuTime(totalTimeInSeconds);
    }
}

OptimizationApplication::~OptimizationApplication()
{
    destroyProblem();

    if(resultFileName)
        delete[] resultFileName;

    destroyHDF5Mutex();

    MPI_Finalize();
}

int OptimizationApplication::getMpiRank()
{
    int mpiRank = -1;

    int mpiInitialized = 0;
    MPI_Initialized(&mpiInitialized);

    if(mpiInitialized) {
        MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    }

    return mpiRank;
}

int OptimizationApplication::getMpiCommSize()
{
    int mpiCommSize = -1;

    int mpiInitialized = 0;
    MPI_Initialized(&mpiInitialized);

    if(mpiInitialized) {
        MPI_Comm_size(MPI_COMM_WORLD, &mpiCommSize);
    }

    return mpiCommSize;
}

