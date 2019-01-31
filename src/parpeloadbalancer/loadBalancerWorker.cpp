#include <parpeloadbalancer/loadBalancerWorker.h>

#ifdef PARPE_ENABLE_MPI

#include <alloca.h>
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <mpi.h>

#define QUEUE_WORKER_H_VERBOSE 0
#define LOADBALANCERWORKER_REPORT_WAITING_TIME 1

#ifdef LOADBALANCERWORKER_REPORT_WAITING_TIME
#include <parpecommon/logging.h>
#endif

namespace parpe {

void LoadBalancerWorker::run(messageHandlerFunc const& messageHandler) {
    bool terminate = false;

    while (!terminate) {
        terminate = waitForAndHandleJobs(messageHandler);
    }
}

bool LoadBalancerWorker::waitForAndHandleJobs(const messageHandlerFunc& messageHandler) {
    int rank, err;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#ifdef LOADBALANCERWORKER_REPORT_WAITING_TIME
    double startTime = MPI_Wtime();
#endif

#if QUEUE_WORKER_H_VERBOSE >= 3
    printf("[%d] Waiting for work.\n", rank);
#endif

    // wait for receiving a single job and check for size
    MPI_Status mpiStatus;
    MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &mpiStatus);
    int msgSize;
    MPI_Get_count(&mpiStatus, MPI_BYTE, &msgSize);
    std::vector<char> buffer(static_cast<unsigned int>(msgSize));

    // receive message
    int source = 0;
    err = MPI_Recv(buffer.data(), msgSize, MPI_BYTE, source, MPI_ANY_TAG,
                   MPI_COMM_WORLD, &mpiStatus);

#if QUEUE_WORKER_H_VERBOSE >= 3
    printf("W%d: Received job %d\n", rank, mpiStatus.MPI_TAG);
#endif
    if (err != MPI_SUCCESS)
        abort();

    if (mpiStatus.MPI_TAG == MPI_TAG_EXIT_SIGNAL) {
        return true;
    }

#ifdef LOADBALANCERWORKER_REPORT_WAITING_TIME
    double endTime = MPI_Wtime();
    double waitedSeconds = (endTime - startTime);
    logmessage(LOGLVL_DEBUG, "Message received after waiting %fs.", rank, waitedSeconds);
#endif

    messageHandler(buffer, mpiStatus.MPI_TAG);

#if QUEUE_WORKER_H_VERBOSE >= 2
    printf("[%d] Job done, sending results, %dB.\n", rank, msgSize);
#endif
    MPI_Send(buffer.data(), buffer.size(), MPI_BYTE, 0, mpiStatus.MPI_TAG, MPI_COMM_WORLD);

    return false;
}

} // namespace parpe

#endif
