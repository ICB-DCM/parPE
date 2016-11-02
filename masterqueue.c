#include "masterqueue.h"
#include <pthread.h>
#include <stdlib.h>
#include <string.h>


typedef struct masterQueueElement_tag {
    queueData *data;
    int jobId;
    // linked list
    struct masterQueueElement_tag *nextElement;
} masterQueueElement;

static masterQueueElement *queueStart = 0;
static masterQueueElement *queueEnd = 0;
static int lastJobId = 0;

static pthread_t queueThread;
static pthread_mutex_t mutexQueue = PTHREAD_MUTEX_INITIALIZER;

// one for each worker, idx 0 not used
static MPI_Request *sendRequests;
static MPI_Request *recvRequests;
static queueData **sentJobs;

static void *masterQueueRun(void *unusedArgument);

static masterQueueElement *getNextJob();

static void sendToWorker(int workerID, masterQueueElement *queueElement);

static void receiveFinished(int workerID);

void initMasterQueue() {
    static bool queueCreated = false;

    // There can only be one queue
    if(!queueCreated) {
        pthread_create(&queueThread, NULL, masterQueueRun, NULL);
        queueCreated = true;
    }
}

// Thread entry point
static void *masterQueueRun(void *unusedArgument) {
    int mpiCommSize;
    MPI_Comm_size(MPI_COMM_WORLD, &mpiCommSize);

    sendRequests = alloca(mpiCommSize * sizeof(MPI_Request));
    recvRequests = alloca(mpiCommSize * sizeof(MPI_Request));
    sentJobs = alloca(mpiCommSize * sizeof(queueData *));

    for(int i = 0; i < mpiCommSize; ++i) // have to init before can wait!
        recvRequests[i] = MPI_REQUEST_NULL;

    while(1) {
        int freeWorkerIndex;
        MPI_Waitany(mpiCommSize, recvRequests, &freeWorkerIndex, MPI_STATUS_IGNORE);

        if(freeWorkerIndex != MPI_UNDEFINED) {
            receiveFinished(freeWorkerIndex);
        } else {
            // no jobs yet, send to first worker
            freeWorkerIndex = 1;
        }

        masterQueueElement *currentQueueElement = getNextJob();

        if(currentQueueElement) {

            sendToWorker(freeWorkerIndex, currentQueueElement);
            sentJobs[freeWorkerIndex] = currentQueueElement->data;
            free(currentQueueElement);
        }

        pthread_yield();
    };

    return 0;
}

static masterQueueElement *getNextJob() {

    pthread_mutex_lock(&mutexQueue);

    masterQueueElement *oldStart = queueStart;

    if(queueStart) {
        queueStart = oldStart->nextElement;
    }

    pthread_mutex_unlock(&mutexQueue);

    return oldStart;
}


static void sendToWorker(int workerID, masterQueueElement *queueElement) {
    int tag = queueElement->jobId;

    queueData *data = queueElement->data;

    MPI_Isend(data->sendBuffer, data->lenSendBuffer, MPI_BYTE, workerID, tag, MPI_COMM_WORLD, &sendRequests[workerID]);

    MPI_Irecv(data->recvBuffer, data->lenRecvBuffer, MPI_BYTE, workerID, tag, MPI_COMM_WORLD, &recvRequests[workerID]);
}


void queueSimulation(queueData *jobData) {
    masterQueueElement *queueElement = malloc(sizeof(masterQueueElement));

    queueElement->data = jobData;
    queueElement->nextElement = 0;

    pthread_mutex_lock(&mutexQueue);

    if(lastJobId == INT_MAX) // Unlikely, but prevent overflow
        lastJobId = 0;

    queueElement->jobId = ++lastJobId;

    if(queueStart) {
        queueEnd->nextElement = queueElement;
    } else {
        queueStart = queueElement;
    }
    queueEnd = queueElement;

    pthread_mutex_unlock(&mutexQueue);
}

void terminateMasterQueue() {
    pthread_mutex_destroy(&mutexQueue);
    pthread_cancel(queueThread);
}


static void receiveFinished(int workerID)
{
    queueData *data = sentJobs[workerID];
    *data->jobDone = true;
}
