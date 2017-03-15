#include "CppUTest/TestHarness_c.h"
#include "CppUTestExt/MockSupport_c.h"

#include <unistd.h>
#include <stdio.h>

#define QUEUE_MASTER_TEST

// include .c because of static functions
#include "queuemaster.c"

TEST_GROUP_C_SETUP(queuemaster) {
    // reset globals
    masterQueue.numWorkers = 0;
    masterQueue.queue = NULL;
    masterQueue.lastJobId = 0;
    masterQueue.mutexQueue = (pthread_mutex_t) PTHREAD_MUTEX_INITIALIZER;
    masterQueue.sendRequests = 0;
    masterQueue.recvRequests = 0;
    masterQueue.sentJobsData = 0;
    memset(&masterQueue.semQueue, 0, sizeof masterQueue.semQueue);
    memset(&masterQueue.queueThread, 0, sizeof masterQueue.queueThread);
}

TEST_GROUP_C_TEARDOWN(queuemaster) {
    mock_c()->checkExpectations();
    mock_c()->clear();
}

// mock
static void assertMPIInitialized() {}

// mock for MPI_Comm_size
int MPI_Comm_size(MPI_Comm comm, int *size) {
    mock_c()->actualCall("MPI_Comm_size");
    *size = 10;
    return MPI_SUCCESS;
}

int MPI_Testany(int count, MPI_Request array_of_requests[], int *index,
                               int *flag, MPI_Status *status) {
    // mock_c()->actualCall("MPI_Testany");

    sleep(1000); // do nothing and wait to be killed

    return 0;
}

TEST_C(queuemaster, test_queueinit) {
    mock_c()->expectOneCall("MPI_Comm_size");
    // Can happen or not, depending on how quick it's terminated
    // mock_c()->expectOneCall("MPI_Testany");

    initMasterQueue();

    CHECK_EQUAL_C_BOOL(true, masterQueue.queue != 0);
    int actVal;
    sem_getvalue(&masterQueue.semQueue, &actVal);
    CHECK_EQUAL_C_BOOL(true, actVal > 0);

    terminateMasterQueue();
}

TEST_C(queuemaster, test_queue) {
    mock_c()->expectNCalls(1, "MPI_Comm_size");
    initMasterQueue();

    JobData data;
    queueSimulation(&data);
    CHECK_EQUAL_C_INT(1, masterQueue.lastJobId);

    JobData data2;
    queueSimulation(&data2);
    CHECK_EQUAL_C_INT(2, masterQueue.lastJobId);
    // getNextJob()

    terminateMasterQueue();
}

TEST_C(queuemaster, test_queue_reinitialization) {
    // should not crash or cause memory leaks
    mock_c()->expectNCalls(1, "MPI_Comm_size");
    initMasterQueue();
    initMasterQueue();
    terminateMasterQueue();
}

TEST_C(queuemaster, test_terminateMasterQueue_noInit) {
    // terminate uninitialized masterQueue should not fail
    terminateMasterQueue();
}
