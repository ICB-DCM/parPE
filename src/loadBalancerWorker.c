#include <string.h>
#include <mpi.h>
#include <assert.h>
#include <alloca.h>
#include <assert.h>
#include <stdbool.h>
#include "loadBalancerWorker.h"
#include "misc.h"

static bool waitForAndHandleJobs(int inMsgSize, int outMesgSize, messageHandlerFp msgHandler, char *buffer);

void loadBalancerWorkerRun(int inMsgSize, int outMesgSize, messageHandlerFp msgHandler) {
    // TODO replace fixed buffer size by MPI_PROBE?
    int bufferSize = (outMesgSize > inMsgSize) ? outMesgSize : inMsgSize;
    char *buffer = alloca(bufferSize);

    bool terminate = false;

    while(!terminate) {
        terminate = waitForAndHandleJobs(inMsgSize, outMesgSize, msgHandler, buffer);
    }

}


bool waitForAndHandleJobs(int inMsgSize, int outMesgSize, messageHandlerFp msgHandler, char *buffer) {
    int rank, err;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

#if QUEUE_WORKER_H_VERBOSE >= 3
    printf("[%d] Waiting for work.\n", rank); fflush(stdout);
#endif

    // wait for receiving a single job
    int source = 0;
    MPI_Status mpiStatus;
    err = MPI_Recv(buffer, inMsgSize, MPI_BYTE, source, MPI_ANY_TAG, MPI_COMM_WORLD, &mpiStatus);

    //printf("W%d: Received job %d\n", rank, mpiStatus.MPI_TAG);
    if(err != MPI_SUCCESS)
        abort();

    if(mpiStatus.MPI_TAG == MPI_TAG_EXIT_SIGNAL)
        return 1;

    msgHandler(buffer, mpiStatus.MPI_TAG);

#if QUEUE_WORKER_H_VERBOSE >= 2
    printf("[%d] Simulation done, sending results (llh: %f). ", rank, result.llh); printDatapath(path); fflush(stdout);
#endif

    MPI_Send(buffer, outMesgSize, MPI_BYTE, 0, mpiStatus.MPI_TAG, MPI_COMM_WORLD);

    return 0;
}
