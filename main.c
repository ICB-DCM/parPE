#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>

#include <mpi.h>
#include <mpe.h>
#include "mpiworker.h"

#include "localoptimization.h"
#include "objectivefunction.h"
#include "logger.h"
#include "masterqueue.h"
#include "masterthread.h"
#include "dataprovider.h"

#include <signal.h>
#include <unistd.h>

volatile sig_atomic_t caughtTerminationSignal = 0;

void getMpeLogIDs();
void describeMpeStates();
void term(int sigNum) { caughtTerminationSignal = 1; }
char *getResultFileName();
void doMasterWork();
void printDebugInfoAndWait();

// MPE event IDs for logging
int mpe_event_begin_simulate, mpe_event_end_simulate;
int mpe_event_begin_getrefs, mpe_event_end_getrefs;
int mpe_event_begin_getdrugs, mpe_event_end_getdrugs;
int mpe_event_begin_aggregate, mpe_event_end_aggregate;

int main(int argc, char **argv)
{
    int mpiCommSize, mpiRank, mpiErr;

    mpiErr = MPI_Init(&argc, &argv);
    mpiErr = MPI_Comm_size(MPI_COMM_WORLD, &mpiCommSize);
    mpiErr = MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    assert(mpiCommSize > 1);

    // printDebugInfoAndWait();

    MPE_Init_log();
    getMpeLogIDs();

    initDataProvider("/home/dweindl/src/CanPathPro-NonMatlab/data/data.h5"); // TODO arguemnt (log: Reading from..

    char *resultFileName = getResultFileName();
    loggerdata datalogger = initResultHDFFile(resultFileName, "/results/multistarts/0"); // TODO
    free(resultFileName);

    if(mpiRank == 0) {

        struct sigaction action;
        action.sa_handler = term;
        sigaction(SIGTERM, &action, NULL);

        describeMpeStates();

        doMasterWork();

        sendTerminationSignalToAllWorkers();

    } else {
        doWorkerWork();
    }

    closeDataProvider();
    closeResultHDFFile(datalogger);

    printf("Finalizing MPE log: mpe.log\n");
    mpiErr = MPE_Finish_log("mpe.log");

    mpiErr = MPI_Finalize();
}

void printDebugInfoAndWait() {
    //int i = 0;
    char hostname[256];
    gethostname(hostname, sizeof(hostname));
    printf("PID %d on %s ready for attach\n", getpid(), hostname);
    fflush(stdout);
    //while (0 == i)
        sleep(15);
}

void getMpeLogIDs() {
    MPE_Log_get_state_eventIDs(&mpe_event_begin_simulate, &mpe_event_end_simulate);
    MPE_Log_get_state_eventIDs(&mpe_event_begin_aggregate, &mpe_event_end_aggregate);
    MPE_Log_get_state_eventIDs(&mpe_event_begin_getrefs, &mpe_event_end_getrefs);
    MPE_Log_get_state_eventIDs(&mpe_event_begin_getdrugs, &mpe_event_end_getdrugs);

}

void describeMpeStates() {
    MPE_Describe_state(mpe_event_begin_simulate, mpe_event_end_simulate, "simulate", "blue:gray");
    MPE_Describe_state(mpe_event_begin_aggregate, mpe_event_end_aggregate, "aggregate", "red:gray");
    MPE_Describe_state(mpe_event_begin_getrefs, mpe_event_end_getrefs, "getrefs", "green:gray");
    MPE_Describe_state(mpe_event_begin_getdrugs, mpe_event_end_getdrugs, "getdrugs", "yellow:gray");
}

char *getResultFileName() {
    time_t timer;
    time(&timer);

    struct tm* tm_info;
    tm_info = localtime(&timer);

    char tmpFileName[50];
    strftime(tmpFileName, 50, "CPP_%Y-%m-%d_%H%M%S_%%05d.h5", tm_info);

    int mpiRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

    char *fileName = malloc(31);
    sprintf(fileName, tmpFileName, mpiRank);

    return fileName;
}

void doMasterWork() {
    initMasterQueue();

    // create threads for multistart batches
    int numMultiStartRuns = getNumMultiStartRuns();
    pthread_t *multiStartThreads = alloca(numMultiStartRuns * sizeof(pthread_t));

    pthread_attr_t threadAttr;
    pthread_attr_init(&threadAttr);
    pthread_attr_setdetachstate(&threadAttr, PTHREAD_CREATE_JOINABLE);

    for(int k = 0; k < numMultiStartRuns; ++k) {
        pthread_create(&multiStartThreads[k], &threadAttr, newMultiStartOptimization, (void *)&k);
    }
    pthread_attr_destroy(&threadAttr);

    // wait for finish
    for(int k = 0; k < numMultiStartRuns; ++k) {
        printf("Waiting for thread k %d\n", k);
        pthread_join(multiStartThreads[k], NULL);
        sleep(0.1);
    }
    printf("All k threads finished.\n");

    terminateMasterQueue();
}
