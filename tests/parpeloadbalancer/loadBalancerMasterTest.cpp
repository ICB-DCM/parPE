#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <parpeloadbalancer/loadBalancerMaster.h>
#include <parpecommon/misc.h>

#include <mpi.h>

#include <cstdio>
#include <unistd.h>

#define QUEUE_MASTER_TEST

using ::testing::_;

static std::function<int(MPI_Comm comm, int *size) > _MPI_Comm_size;

// mock for MPI_Comm_size
int MPI_Comm_size(MPI_Comm comm, int *size) {
    _MPI_Comm_size(comm, size);
    *size = 10;
    return MPI_SUCCESS;
}

int MPI_Testany(int count, MPI_Request array_of_requests[], int *index,
                int *flag, MPI_Status *status) {
    sleep(1000); // do nothing and wait to be killed

    return 0;
}

int MPI_Iprobe(int source, int tag, MPI_Comm comm, int *flag,
               MPI_Status *status) {
    sleep(1000);
    return 0;
}


class MockMPI {
public:
    MOCK_CONST_METHOD2(MPI_Comm_size, int(MPI_Comm comm, int *size));
    MOCK_CONST_METHOD5(MPI_Testany, int(int count,
                                        MPI_Request array_of_requests[],
                                        int *index,
                                        int *flag, MPI_Status *status));
    MOCK_CONST_METHOD5(MPI_Iprobe, int(int source, int tag, MPI_Comm comm,
                                       int *flag, MPI_Status *status));

    MockMPI() {
        _MPI_Comm_size = [this](MPI_Comm comm, int *size){ return MPI_Comm_size(comm, size); };
    }
};

class queuemaster : public ::testing::Test {

protected:
    MockMPI mockMpi;
};


#include <loadBalancerMaster.cpp>

TEST_F(queuemaster, test_queueinit) {
    EXPECT_CALL(mockMpi, MPI_Comm_size(_, _)).Times(1);
    // Can happen or not, depending on how quick it's terminated
    // mock().expectOneCall("MPI_Testany");

    parpe::LoadBalancerMaster lbm;
    lbm.run();
    EXPECT_EQ(lbm.getNumQueuedJobs(), 0);
    lbm.terminate();
}

TEST_F(queuemaster, test_queue) {
    EXPECT_CALL(mockMpi, MPI_Comm_size(_, _)).Times(1);
    parpe::LoadBalancerMaster lbm;
    lbm.run();

    EXPECT_EQ(lbm.getNumQueuedJobs(), 0);

    parpe::JobData data;
    lbm.queueJob(&data);
    EXPECT_EQ(lbm.getNumQueuedJobs(), 1);

    parpe::JobData data2;
    lbm.queueJob(&data2);
    EXPECT_EQ(lbm.getNumQueuedJobs(), 2);

    lbm.terminate();
}

TEST_F(queuemaster, test_terminateMasterQueue_noInit) {
    // terminate uninitialized masterQueue should not fail
    parpe::LoadBalancerMaster lbm;
    lbm.terminate();
}
