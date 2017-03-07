#include "CppUTest/TestHarness_c.h"
#include "CppUTestExt/MockSupport_c.h"

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
}

// TODO: requires MPI_Init
//TEST_C(masterqueue, test_queueinit) {
//    initMasterQueue();

//    CHECK_EQUAL_C_BOOL(true, queueCreated);
//    int actVal;
//    sem_getvalue(&actValsemQueue, &actVal);
//    CHECK_EQUAL_C_BOOL(true, actVal > 0);
//}

//TEST_C(masterqueue, test_queue) {

//}

TEST_C(masterqueue, test_terminateMasterQueue_noInit) {
    // terminate uninitialized masterQueue should not fail
    terminateMasterQueue();
}
