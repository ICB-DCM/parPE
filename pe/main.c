#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <pthread.h>
#include <mpi.h>
#include <getopt.h>

#include "loadBalancerMaster.h"
#include "simulationWorker.h"
#include "resultwriter.h"
#include "masterthread.h"
#include "dataprovider.h"
#include "misc.h"

#undef INSTALL_SIGNAL_HANDLER
#ifdef INSTALL_SIGNAL_HANDLER
#include <signal.h>
volatile sig_atomic_t caughtTerminationSignal = 0;
void term(int sigNum) { caughtTerminationSignal = 1; }
#endif

#ifdef USE_MPE
#include <mpe.h>
// MPE event IDs for logging
int mpe_event_begin_simulate, mpe_event_end_simulate;
int mpe_event_begin_getrefs, mpe_event_end_getrefs;
int mpe_event_begin_getdrugs, mpe_event_end_getdrugs;
int mpe_event_begin_aggregate, mpe_event_end_aggregate;
#endif

// global mutex for HDF5 library calls
pthread_mutex_t mutexHDF;

typedef enum operationType_tag {OP_TYPE_PARAMETER_ESTIMATION, OP_TYPE_GRADIENT_CHECK} operationTypeEnum;

operationTypeEnum opType = OP_TYPE_PARAMETER_ESTIMATION;

optimizerEnum optimizer = OPTIMIZER_IPOPT;

// program options
const char *shortOptions = "dhvt:o:";
static struct option const longOptions[] = {
    {"debug", no_argument, NULL, 'd'},
    {"print-worklist", no_argument, NULL, 'p'},
    {"help", no_argument, NULL, 'h'},
    {"version", no_argument, NULL, 'v'},
    {"task", required_argument, NULL, 't'},
    {"optimizer", required_argument, NULL, 'o'},
    {NULL, 0, NULL, 0}
};

// HDF5 file for reading options and data
static const char *inputFile;

// intercept malloc calls & check results
extern void *__libc_malloc(size_t);

inline void *malloc(size_t size)
{
    void *p = __libc_malloc(size);

    if(size && !p) {
        fprintf(stderr, "Could not allocate memory.\n");
        abort();
    }

    return p;
}

void init(int argc, char **argv);

void doMasterWork();

int finalize(clock_t begin);

void finalizeTiming(clock_t begin);

#ifdef USE_MPE
void getMpeLogIDs();
#endif

void describeMpeStates();

char *getResultFileName();

void initHDF5Mutex();

int parseOptions(int argc, char **argv);

void sendTerminationSignalToAllWorkers();

int main(int argc, char **argv)
{
    clock_t begin = clock();

    init(argc, argv);

    int mpiRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

    if(mpiRank == 0) {
        doMasterWork();
    } else {
        doWorkerWork();
    }

    return finalize(begin);
}

void init(int argc, char **argv) {
    // Seed random number generator
    srand(1337);

    int status;
    status = parseOptions(argc, argv);
    if(status)
        exit(status);

    int mpiErr = MPI_Init(&argc, &argv);
    if(mpiErr != MPI_SUCCESS) {
        logmessage(LOGLVL_CRITICAL, "Problem initializing MPI. Exiting.");
        exit(1);
    }

    int mpiCommSize;
    MPI_Comm_size(MPI_COMM_WORLD, &mpiCommSize);

    if(mpiCommSize < 2) {
        logmessage(LOGLVL_CRITICAL, "Need at least 2 MPI processes, but got %d %d. Exiting.", mpiCommSize);
        exit(1);
    }

    // printDebugInfoAndWait();
    printMPIInfo();

#ifdef USE_MPE
    MPE_Init_log();
    getMpeLogIDs();
#endif

    initHDF5Mutex();

    logmessage(LOGLVL_INFO, "Reading options and data from '%s'.", inputFile);
    status = initDataProvider(inputFile);

    if(status != 0)
        exit(status);

    char *resultFileName = getResultFileName();
    status = initResultHDFFile(resultFileName);
    free(resultFileName);

    if(status != 0)
        exit(status);
}


void doMasterWork() {

    dataproviderPrintInfo();

#ifdef INSTALL_SIGNAL_HANDLER
    struct sigaction action;
    action.sa_handler = term;
    sigaction(SIGTERM, &action, NULL);
#endif

#ifdef USE_MPE
    describeMpeStates();
#endif

    switch (opType) {
    case OP_TYPE_GRADIENT_CHECK:
        startObjectiveFunctionGradientCheck();
        break;
    default:
        startParameterEstimation(optimizer);
        break;
    }

    sendTerminationSignalToAllWorkers();
}


int finalize(clock_t begin) {
    finalizeTiming(begin);

    closeDataProvider();

    closeResultHDFFile();

    pthread_mutex_destroy(&mutexHDF);

    logProcessStats();

#ifdef USE_MPE
    logmessage(LOGLVL_DEBUG, "Finalizing MPE log: mpe.log");
    MPE_Finish_log("mpe.log");
#endif

    logmessage(LOGLVL_DEBUG, "Finalizing MPI");
    MPI_Finalize();

    return 0;
}


void finalizeTiming(clock_t begin) {
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
        saveTotalWalltime(totalTimeInSeconds);
    }
}

#ifdef USE_MPE
void getMpeLogIDs() {
    MPE_Log_get_state_eventIDs(&mpe_event_begin_simulate, &mpe_event_end_simulate);
    MPE_Log_get_state_eventIDs(&mpe_event_begin_aggregate, &mpe_event_end_aggregate);
    MPE_Log_get_state_eventIDs(&mpe_event_begin_getrefs, &mpe_event_end_getrefs);
    MPE_Log_get_state_eventIDs(&mpe_event_begin_getdrugs, &mpe_event_end_getdrugs);
}
#endif


#ifdef USE_MPE
void describeMpeStates() {
    MPE_Describe_state(mpe_event_begin_simulate,  mpe_event_end_simulate,  "simulate",  "blue:gray");
    MPE_Describe_state(mpe_event_begin_aggregate, mpe_event_end_aggregate, "aggregate", "red:gray");
    MPE_Describe_state(mpe_event_begin_getrefs,   mpe_event_end_getrefs,   "getrefs",   "green:gray");
    MPE_Describe_state(mpe_event_begin_getdrugs,  mpe_event_end_getdrugs,  "getdrugs",  "yellow:gray");
}
#endif

char *getResultFileName() {
    // create directory for each compute node
    char procName[MPI_MAX_PROCESSOR_NAME];
    int procNameLen;
    MPI_Get_processor_name(procName, &procNameLen);

    createDirectoryIfNotExists(procName);

    char tmpFileName[200];
    strFormatCurrentLocaltime(tmpFileName, 200, "%%s/CPP_%Y-%m-%d_%H%M%S_%%05d.h5");

    int mpiRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

    char *fileName = malloc(MPI_MAX_PROCESSOR_NAME + 32);
    sprintf(fileName, tmpFileName, procName, mpiRank);

    return fileName;
}


void initHDF5Mutex() {
    pthread_mutexattr_t attr;
    pthread_mutexattr_init(&attr);
    pthread_mutexattr_settype(&attr, PTHREAD_MUTEX_RECURSIVE);
    pthread_mutex_init(&mutexHDF, &attr);
    pthread_mutexattr_destroy(&attr);

    H5dont_atexit();
}

int parseOptions (int argc, char **argv) {

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
            if(strcmp(optarg, "ceres") == 0)
                optimizer = OPTIMIZER_CERES;
            break;
        case 'v':
            printf("Version: %s\n", GIT_VERSION);
            return 1;
        case 'h':
            printf("Usage: %s [OPTION]... FILE\n", argv[0]);
            printf("FILE: HDF5 data file");
            printf("Options: \n"
                   "  -o, --optimizer Which optimizer to use Ipopt (default) Ceres (not yet implemented)"
                   "  -t, --task    What to do? Parameter estimation (default) or check gradient ('gradient_check')"
                   "  -h, --help    Print this help text\n"
                   "  -v, --version Print version info\n"
                   );
            return 1;
        default:
            printf("Unrecognized option: %c\n", c);
        }
    }

    if(optind < argc) {
        inputFile = argv[optind++];
    } else {
        logmessage(LOGLVL_CRITICAL, "Must provide input file as first and only argument to %s.", argv[0]);
        return 1;
    }

    return 0;
}


void sendTerminationSignalToAllWorkers()
{
    int commSize;
    MPI_Comm_size(MPI_COMM_WORLD, &commSize);

    MPI_Request reqs[commSize - 1];

    for(int i = 1; i < commSize; ++i) {
        reqs[i - 1] =  MPI_REQUEST_NULL;
        MPI_Isend(MPI_BOTTOM, 0, MPI_INT, i, 0, MPI_COMM_WORLD, &reqs[i - 1]);
    }
    logmessage(LOGLVL_INFO, "Sent termination signal to workers.");
    MPI_Waitall(commSize - 1, reqs, MPI_STATUS_IGNORE);
}
