#include <bits/stl_tree.h>

#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"

#include <cstdio>
#include <mpi.h>
#include <unistd.h>
#include <misc.h>

#define QUEUE_MASTER_TEST

TEST_GROUP(queuemaster){

    void setup(){

    }

    void teardown(){mock().checkExpectations();
mock().clear();
}
}
;

// mock
static void assertMPIInitialized() {}

// mock for MPI_Comm_size
int MPI_Comm_size(MPI_Comm comm, int *size) {
    mock().actualCall("MPI_Comm_size");
    *size = 10;
    return MPI_SUCCESS;
}

int MPI_Testany(int count, MPI_Request array_of_requests[], int *index,
                int *flag, MPI_Status *status) {
    // mock().actualCall("MPI_Testany");

    sleep(1000); // do nothing and wait to be killed

    return 0;
}

int MPI_Iprobe(int source, int tag, MPI_Comm comm, int *flag,
               MPI_Status *status) {
    sleep(1000);
    return 0;
}

#include <LoadBalancerMaster.cpp>

TEST(queuemaster, test_queueinit) {
    mock().expectOneCall("MPI_Comm_size");
    // Can happen or not, depending on how quick it's terminated
    // mock().expectOneCall("MPI_Testany");

    parpe::LoadBalancerMaster lbm;
    lbm.run();

    //    CHECK_C(loadBalancer.queue != 0);
    //    int actVal;
    //    sem_getvalue(&loadBalancer.semQueue, &actVal);
    //    CHECK_C(actVal > 0);

    lbm.terminate();
}

TEST(queuemaster, test_queue) {
    mock().expectNCalls(1, "MPI_Comm_size");
    parpe::LoadBalancerMaster lbm;
    lbm.run();

    parpe::JobData data;
    lbm.queueJob(&data);
    //    CHECK_EQUAL(1, lbm.lastJobId);

    parpe::JobData data2;
    lbm.queueJob(&data2);
    //    CHECK_EQUAL(2, loadBalancer.lastJobId);
    // getNextJob()

    lbm.terminate();
}

//
// TEST(queuemaster, test_terminateMasterQueue_noInit) {
//    // terminate uninitialized masterQueue should not fail
//    loadBalancerTerminate();
//}
