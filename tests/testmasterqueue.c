#include "CppUTest/TestHarness_c.h"
#include "CppUTestExt/MockSupport_c.h"

#include <unistd.h>

// include .c because of static functions
#include "masterqueue.c"

TEST_GROUP_C_SETUP(masterqueue) {
    // reset globals
    queueCreated = false;
    queueStart = 0;
    queueEnd = 0;
    lastJobId = 0;
    memset(&semQueue, 0, sizeof semQueue);
    memset(&queueThread, 0, sizeof queueThread);
    mutexQueue = (pthread_mutex_t) PTHREAD_MUTEX_INITIALIZER;
    sendRequests = 0;
    recvRequests = 0;
    sentJobsData = 0;
}

TEST_GROUP_C_TEARDOWN(masterqueue) {
    mock_c()->checkExpectations();
    mock_c()->clear();
}

// mock for MPI_Comm_size
int MPI_Comm_size(MPI_Comm comm, int *size) {
    mock_c()->actualCall("MPI_Comm_size");
    *size = 10;
    return 1; // ?
}

int MPI_Testany(int count, MPI_Request array_of_requests[], int *index,
                               int *flag, MPI_Status *status) {
    // mock_c()->actualCall("MPI_Testany");

    sleep(1000); // do nothing and wait to be killed

    return 0;
}

TEST_C(masterqueue, test_queueinit) {
    mock_c()->expectOneCall("MPI_Comm_size");
    // Can happen or not, depending on how quick it's terminated
    // mock_c()->expectOneCall("MPI_Testany");

    initMasterQueue();

    CHECK_EQUAL_C_BOOL(true, queueCreated);
    int actVal;
    sem_getvalue(&semQueue, &actVal);
    CHECK_EQUAL_C_BOOL(true, actVal > 0);

    terminateMasterQueue();
}

TEST_C(masterqueue, test_queue) {
    mock_c()->expectNCalls(1, "MPI_Comm_size");
    initMasterQueue();

    queueData data;
    queueSimulation(&data);
    CHECK_EQUAL_C_BOOL(false, queueEnd == 0);
    CHECK_EQUAL_C_BOOL(false, queueStart == 0);
    CHECK_EQUAL_C_BOOL(true, queueEnd == queueStart);
    CHECK_EQUAL_C_INT(1, lastJobId);

    queueSimulation(&data);
    CHECK_EQUAL_C_BOOL(false, queueEnd == queueStart);
    CHECK_EQUAL_C_INT(2, lastJobId);
    // getNextJob()

    terminateMasterQueue();
}

TEST_C(masterqueue, test_queue_reinitialization) {
    // should not crash or cause memory leaks
    mock_c()->expectNCalls(1, "MPI_Comm_size");
    initMasterQueue();
    initMasterQueue();
    terminateMasterQueue();
}


TEST_C(masterqueue, test_terminateMasterQueue_noInit) {
    // terminate uninitialized masterQueue should not fail
    terminateMasterQueue();
}
