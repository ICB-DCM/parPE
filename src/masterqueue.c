#include "masterqueue.h"
#include "misc.h"

#include <pthread.h>
#include <semaphore.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdbool.h>
#include <limits.h>

typedef struct masterQueueElement_tag {
    queueData *data;
    int jobId;
    // linked list
    struct masterQueueElement_tag *nextElement;
} masterQueueElement;

static int numWorkers;
static bool queueCreated = false;
static masterQueueElement *queueStart = 0;
static masterQueueElement *queueEnd = 0;
static int lastJobId = 0;

static sem_t semQueue;
static pthread_t queueThread;
static pthread_mutex_t mutexQueue = PTHREAD_MUTEX_INITIALIZER;

// one for each worker, index is off by one from MPI rank
static MPI_Request *sendRequests = 0;
static MPI_Request *recvRequests = 0;
static queueData **sentJobsData = 0;

static void *masterQueueRun(void *unusedArguemnt);

static masterQueueElement *getNextJob();

static void sendToWorker(int workerIdx, masterQueueElement *queueElement);

static void receiveFinished(int workerID, int jobID);

static void freeQueueElements();

void initMasterQueue() {
    // There can only be one queue
    if(!queueCreated) {
        int mpiCommSize;
        MPI_Comm_size(MPI_COMM_WORLD, &mpiCommSize);

        numWorkers = mpiCommSize - 1;
        sendRequests = malloc(numWorkers * sizeof(MPI_Request));
        recvRequests = malloc(numWorkers * sizeof(MPI_Request));
        sentJobsData = malloc(numWorkers * sizeof(queueData *));

        for(int i = 0; i < numWorkers; ++i) // have to initialize before can wait!
            recvRequests[i] = MPI_REQUEST_NULL;

        // Create semaphore to limit queue length
        // and avoid huge memory allocation for all send and receive buffers
        unsigned int queueMaxLength = mpiCommSize;
#ifdef SEM_VALUE_MAX
        queueMaxLength = SEM_VALUE_MAX < queueMaxLength ? SEM_VALUE_MAX : queueMaxLength;
#endif
        sem_init(&semQueue, 0, queueMaxLength);
        pthread_create(&queueThread, NULL, masterQueueRun, 0);
        queueCreated = true;
    }
}

// Thread entry point
static void *masterQueueRun(void *unusedArguemnt) {

    // dispatch queued work packages
    while(1) {
        // check if any job finished
        MPI_Status status;
        int finishedWorkerIdx = 0;

        // handle all finished jobs, if any
        while(1) {
            int dummy;
            MPI_Testany(numWorkers, recvRequests, &finishedWorkerIdx, &dummy, &status);

            if(finishedWorkerIdx >= 0) {
                // dummy == 1 despite finishedWorkerIdx == MPI_UNDEFINED
                // some job is finished
                receiveFinished(finishedWorkerIdx, status.MPI_TAG);
            } else {
                // there was nothing to be finished
                break;
            }
        }

        // getNextFreeWorker
        int freeWorkerIndex = finishedWorkerIdx;

        if(freeWorkerIndex < 0) { // no job finished recently, checked free slots
            for(int i = 0; i < numWorkers; ++i) {
                if(recvRequests[i] == MPI_REQUEST_NULL) {
                    freeWorkerIndex = i;
                    break;
                }
            }
        }

        if(freeWorkerIndex < 0) {
            // all workers are busy, wait for next one to finish
            MPI_Status status;
            MPI_Waitany(numWorkers, recvRequests, &freeWorkerIndex, &status);

            assert(freeWorkerIndex != MPI_UNDEFINED);
            receiveFinished(freeWorkerIndex, status.MPI_TAG);
        }

        masterQueueElement *currentQueueElement = getNextJob();

        if(currentQueueElement) {
            sendToWorker(freeWorkerIndex, currentQueueElement);
            sentJobsData[freeWorkerIndex] = currentQueueElement->data;
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


static void sendToWorker(int workerIdx, masterQueueElement *queueElement) {
    int tag = queueElement->jobId;
    int workerRank = workerIdx + 1;
    queueData *data = queueElement->data;

#ifdef MASTER_QUEUE_H_SHOW_COMMUNICATION
    printf("\x1b[36mSent job %d to %d.\x1b[0m\n", tag, workerRank);
#endif

    MPI_Isend(data->sendBuffer, data->lenSendBuffer, MPI_BYTE, workerRank, tag, MPI_COMM_WORLD, &sendRequests[workerIdx]);

    MPI_Irecv(data->recvBuffer, data->lenRecvBuffer, MPI_BYTE, workerRank, tag, MPI_COMM_WORLD, &recvRequests[workerIdx]);

    sem_post(&semQueue);
}


void queueSimulation(queueData *jobData) {
    assert(queueCreated);

    sem_wait(&semQueue);

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
    sem_destroy(&semQueue);

    if(sentJobsData)
        free(sentJobsData);
    if(sendRequests)
        free(sendRequests);
    if(recvRequests)
        free(recvRequests);

    freeQueueElements();
}

/**
 * @brief freeQueueElements Deallocate queueElement memory when non-empty queue is terminated.
 * Freeing masterQueueElement.queueData memory is the users problem.
 */
static void freeQueueElements() {
    masterQueueElement *nextElement = queueStart;
    while(nextElement) {
        masterQueueElement *curElement = nextElement;
        nextElement = curElement->nextElement;
        free(curElement);
    }
}

static void receiveFinished(int workerID, int jobID) {

#ifdef MASTER_QUEUE_H_SHOW_COMMUNICATION
    printf("\x1b[32mReceived result for job %d from %d\x1b[0m\n", jobID, workerID + 1);
#endif

    queueData *data = sentJobsData[workerID];
    ++(*data->jobDone);
    pthread_mutex_lock(data->jobDoneChangedMutex);
    pthread_cond_signal(data->jobDoneChangedCondition);
    pthread_mutex_unlock(data->jobDoneChangedMutex);

    recvRequests[workerID] = MPI_REQUEST_NULL;
}
