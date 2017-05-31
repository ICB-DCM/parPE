#include "CppUTest/TestHarness_c.h"
#include "CppUTestExt/MockSupport_c.h"

#include <unistd.h>
#include <stdio.h>

#define QUEUE_MASTER_TEST

// include .c because of static functions
#include "loadBalancerMaster.c"

TEST_GROUP_C_SETUP(queuemaster) {
    // reset globals
    loadBalancer.numWorkers = 0;
    loadBalancer.queue = NULL;
    loadBalancer.lastJobId = 0;
    loadBalancer.mutexQueue = (pthread_mutex_t) PTHREAD_MUTEX_INITIALIZER;
    loadBalancer.sendRequests = 0;
    loadBalancer.workerIsBusy = 0;
    loadBalancer.sentJobsData = 0;
    memset(&loadBalancer.semQueue, 0, sizeof loadBalancer.semQueue);
    memset(&loadBalancer.queueThread, 0, sizeof loadBalancer.queueThread);
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

    loadBalancerStartMaster();

    CHECK_C(loadBalancer.queue != 0);
    int actVal;
    sem_getvalue(&loadBalancer.semQueue, &actVal);
    CHECK_C(actVal > 0);

    loadBalancerTerminate();
}

TEST_C(queuemaster, test_queue) {
    mock_c()->expectNCalls(1, "MPI_Comm_size");
    loadBalancerStartMaster();

    JobData data;
    loadBalancerQueueJob(&data);
    CHECK_EQUAL_C_INT(1, loadBalancer.lastJobId);

    JobData data2;
    loadBalancerQueueJob(&data2);
    CHECK_EQUAL_C_INT(2, loadBalancer.lastJobId);
    // getNextJob()

    loadBalancerTerminate();
}

TEST_C(queuemaster, test_queue_reinitialization) {
    // should not crash or cause memory leaks
    mock_c()->expectNCalls(1, "MPI_Comm_size");
    loadBalancerStartMaster();
    loadBalancerStartMaster();
    loadBalancerTerminate();
}

TEST_C(queuemaster, test_terminateMasterQueue_noInit) {
    // terminate uninitialized masterQueue should not fail
    loadBalancerTerminate();
}
