#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <alloca.h>
#include <pthread.h>
#include <mpi.h>
#include <mpe.h>
#include <getopt.h>

#undef INSTALL_SIGNAL_HANDLER
#ifdef INSTALL_SIGNAL_HANDLER
#include <signal.h>
#endif

#include <include/amici.h>

#include "mpiworker.h"
#include "localoptimization.h"
#include "objectivefunction.h"
#include "resultwriter.h"
#include "masterqueue.h"
#include "masterthread.h"
#include "dataprovider.h"
#include "misc.h"

#ifdef INSTALL_SIGNAL_HANDLER
volatile sig_atomic_t caughtTerminationSignal = 0;
void term(int sigNum) { caughtTerminationSignal = 1; }
#endif

void printMPIInfo();
void getMpeLogIDs();
void describeMpeStates();
char *getResultFileName();
void startParameterEstimation();
void printDebugInfoAndWait();
void initHDF5Mutex();
int parseOptions(int argc, char **argv);

// MPE event IDs for logging
int mpe_event_begin_simulate, mpe_event_end_simulate;
int mpe_event_begin_getrefs, mpe_event_end_getrefs;
int mpe_event_begin_getdrugs, mpe_event_end_getdrugs;
int mpe_event_begin_aggregate, mpe_event_end_aggregate;

// global mutex for HDF5 library calls
pthread_mutex_t mutexHDF;

static struct option const long_options[] = {
    {"debug", no_argument, NULL, 'd'},
    {"print-worklist", no_argument, NULL, 'p'},
    {"help", no_argument, NULL, 'h'},
    {"version", no_argument, NULL, 'v'},
    {NULL, 0, NULL, 0}
};

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

int benchmarkObj(const double *theta, double *llh) {
    int status;
    datapath dp;
    dp.idxExperiment = 1;
    dp.idxGenotype = 1;
    dp.idxLocalOptimization = 0;
    dp.idxLocalOptimizationIteration = 0;
    dp.idxMultiStart = 0;
    UserData *udata = getMyUserData();
    //int iterationsDone = 0;

    udata->am_sensi_meth = 0;
    //ReturnData *rdata = getSteadystateSolutionForExperiment(dp, udata, &status, 0, &iterationsDone);

    memcpy(udata->am_p, theta, sizeof(double) * NUM_OPTIMIZATION_PARAMS);
    ExpData *edata = getExperimentalDataForExperiment(dp, udata);

    ReturnData *rdata = getSimulationResults(udata, edata, &status);

    *llh = *rdata->am_llhdata;

    myFreeExpData(edata);
    freeUserDataC(udata);
    freeReturnData(rdata);

    return 0;
}


int benchmarkObjGrad(const double *theta, double *gradient) {
    int status;
    datapath dp;
    dp.idxExperiment = 1;
    dp.idxGenotype = 1;
    dp.idxLocalOptimization = 0;
    dp.idxLocalOptimizationIteration = 0;
    dp.idxMultiStart = 0;
    UserData *udata = getMyUserData();
    udata->am_sensi_meth = 2;

    //int iterationsDone = 0;
    //ReturnData *rdata = getSteadystateSolutionForExperiment(dp, udata, &status, 0, &iterationsDone);

    memcpy(udata->am_p, theta, sizeof(double) * NUM_OPTIMIZATION_PARAMS);
    ExpData *edata = getExperimentalDataForExperiment(dp, udata);

    ReturnData *rdata = getSimulationResults(udata, edata, &status);

    memcpy(gradient, rdata->am_sllhdata, sizeof(double) * NUM_OPTIMIZATION_PARAMS);

    myFreeExpData(edata);
    freeUserDataC(udata);
    freeReturnData(rdata);

    return 0;
}

int main(int argc, char **argv)
{
    int status = 0;

    status = parseOptions(argc, argv);
    if(status)
        return status;

    initHDF5Mutex();
    logmessage(LOGLVL_INFO, "Reading options and data from '%s'.", inputFile);
    status = initDataProvider(inputFile);
    datapath dp;
    dp.idxExperiment = 1;
    dp.idxGenotype = 1;
    dp.idxLocalOptimization = 0;
    dp.idxLocalOptimizationIteration = 0;
    dp.idxMultiStart = 0;

    double theta[NUM_OPTIMIZATION_PARAMS];
    getRandomInitialThetaFromFile(dp, theta, AMI_SCALING_LOG10);
    int indices[NUM_OPTIMIZATION_PARAMS];
    for(int i = 0; i < NUM_OPTIMIZATION_PARAMS; ++i)
        indices[i] = i;
    for(int i = 0; i < NUM_OPTIMIZATION_PARAMS; ++i) {
        int a = rand() / (double) RAND_MAX * NUM_OPTIMIZATION_PARAMS;
        int tmp = indices[a];
        indices[a] = indices[i];
        indices[i] = tmp;
    }

    checkGradient(benchmarkObj, benchmarkObjGrad, NUM_OPTIMIZATION_PARAMS, theta, 1e-8, indices, NUM_OPTIMIZATION_PARAMS);

    exit(0);

    // Seed random number generator
    srand(1337);

    int mpiCommSize, mpiRank, mpiErr;

    status = parseOptions(argc, argv);
    if(status)
        return status;

    mpiErr = MPI_Init(&argc, &argv);
    if(mpiErr != MPI_SUCCESS) {
        logmessage(LOGLVL_CRITICAL, "Problem initializing MPI. Exiting.");
        exit(1);
    }

    MPI_Comm_size(MPI_COMM_WORLD, &mpiCommSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    if(mpiCommSize < 2) {
        logmessage(LOGLVL_CRITICAL, "Need at least 2 MPI processes. Exiting.");
        exit(1);
    }

    double startTime = MPI_Wtime();

    // printDebugInfoAndWait();
    printMPIInfo();

    MPE_Init_log();
    getMpeLogIDs();

    initHDF5Mutex();

    logmessage(LOGLVL_INFO, "Reading options and data from '%s'.", inputFile);
    status = initDataProvider(inputFile);

    if(status != 0)
        exit(1);

    char *resultFileName = getResultFileName();
    status = initResultHDFFile(resultFileName);
    free(resultFileName);
    if(status != 0)
        exit(1);

    if(mpiRank == 0) {
        dataproviderPrintInfo();

#ifdef INSTALL_SIGNAL_HANDLER
        struct sigaction action;
        action.sa_handler = term;
        sigaction(SIGTERM, &action, NULL);
#endif
        describeMpeStates();

        startParameterEstimation();

        sendTerminationSignalToAllWorkers();

    } else {
        doWorkerWork();
    }

    closeDataProvider();

    double endTime = MPI_Wtime();
    double wallTimeSeconds = (endTime - startTime);

    double totalTimeInSeconds = 0;
    MPI_Reduce(&wallTimeSeconds, &totalTimeInSeconds, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if(mpiRank == 0) {
        logmessage(LOGLVL_INFO, "Walltime: %fs, total compute time:%fs", wallTimeSeconds, totalTimeInSeconds);
        saveTotalWalltime(totalTimeInSeconds);
    }

    closeResultHDFFile();

    pthread_mutex_destroy(&mutexHDF);

    logProcessStats();

    logmessage(LOGLVL_DEBUG, "Finalizing MPE log: mpe.log");
    MPE_Finish_log("mpe.log");

    logmessage(LOGLVL_DEBUG, "Finalizing MPI");
    MPI_Finalize();
}

void getMpeLogIDs() {
    MPE_Log_get_state_eventIDs(&mpe_event_begin_simulate, &mpe_event_end_simulate);
    MPE_Log_get_state_eventIDs(&mpe_event_begin_aggregate, &mpe_event_end_aggregate);
    MPE_Log_get_state_eventIDs(&mpe_event_begin_getrefs, &mpe_event_end_getrefs);
    MPE_Log_get_state_eventIDs(&mpe_event_begin_getdrugs, &mpe_event_end_getdrugs);
}

void describeMpeStates() {
    MPE_Describe_state(mpe_event_begin_simulate,  mpe_event_end_simulate,  "simulate",  "blue:gray");
    MPE_Describe_state(mpe_event_begin_aggregate, mpe_event_end_aggregate, "aggregate", "red:gray");
    MPE_Describe_state(mpe_event_begin_getrefs,   mpe_event_end_getrefs,   "getrefs",   "green:gray");
    MPE_Describe_state(mpe_event_begin_getdrugs,  mpe_event_end_getdrugs,  "getdrugs",  "yellow:gray");
}

char *getResultFileName() {
    // create directory for each compute node
    char procName[MPI_MAX_PROCESSOR_NAME];
    int procNameLen;
    MPI_Get_processor_name(procName, &procNameLen);

    struct stat st = {0};

    if (stat(procName, &st) == -1) {
        mkdir(procName, 0700);
    }

    // generate file name
    time_t timer;
    time(&timer);

    struct tm* tm_info;
    tm_info = localtime(&timer);

    char tmpFileName[200];
    strftime(tmpFileName, 200, "%%s/CPP_%Y-%m-%d_%H%M%S_%%05d.h5", tm_info);

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
        int option_index = 0;
        c = getopt_long (argc, argv, "dvh", long_options, &option_index);

        if (c == -1)
            break;

        switch (c) {
        case 'd':
            printDebugInfoAndWait();
        case 'v':
            printf("Version: %s\n", GIT_VERSION);
            break;
        case 'h':
            printf("Usage: %s [OPTION]... FILE\n", argv[0]);
            printf("FILE: HDF5 data file");
            printf("Options: \n"
                   "  -h, --help    Print this help text\n"
                   "  -v, --version Print version info\n"
                   );
            break;
        default:
            printf("%c\n", c);
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
