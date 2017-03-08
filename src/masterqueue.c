#include "masterqueue.h"
#include "misc.h"

#include <pthread.h>
#include <semaphore.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdbool.h>
#include <limits.h>

/** a new object is created upon queueSimulation() */
typedef struct masterQueueElement_tag {
    queueData *data;
    int jobId;
    // linked list
    struct masterQueueElement_tag *nextElement;
} masterQueueElement;

typedef struct masterQueue_tag {
    int numWorkers;
    bool queueCreated;
    masterQueueElement *queueStart;
    masterQueueElement *queueEnd;
    int lastJobId;
    // one for each worker, index is off by one from MPI rank
    // because no job is sent to master (rank 0)
    MPI_Request *sendRequests;
    MPI_Request *recvRequests;
    queueData **sentJobsData;
    pthread_mutex_t mutexQueue;
    sem_t semQueue;
    pthread_t queueThread;
} masterQueueStruct;

#define MASTERQUEUE_QUEUE_INITIALIZER {0, false, NULL, NULL, 0, NULL, NULL, NULL, PTHREAD_MUTEX_INITIALIZER, {}, 0}


static masterQueueStruct masterQueue = MASTERQUEUE_QUEUE_INITIALIZER;


static void *masterQueueRun(void *unusedArgument);

static void assertMPIInitialized();

static masterQueueElement *getNextJob();

static void sendToWorker(int workerIdx, masterQueueElement *queueElement);

static void masterQueueAppendElement(masterQueueElement *queueElement);

static void receiveFinished(int workerID, int jobID);

static void freeQueueElements();


/**
 * @brief initMasterQueue Intialize queue.
 */
void initMasterQueue() {
    // There can only be one queue
    if(!masterQueue.queueCreated) {
        assertMPIInitialized();

        int mpiCommSize;
        MPI_Comm_size(MPI_COMM_WORLD, &mpiCommSize);

        masterQueue.numWorkers = mpiCommSize - 1;
        masterQueue.sendRequests = malloc(masterQueue.numWorkers * sizeof(MPI_Request));
        masterQueue.recvRequests = malloc(masterQueue.numWorkers * sizeof(MPI_Request));
        masterQueue.sentJobsData = malloc(masterQueue.numWorkers * sizeof(queueData *));

        for(int i = 0; i < masterQueue.numWorkers; ++i) // have to initialize before can wait!
            masterQueue.recvRequests[i] = MPI_REQUEST_NULL;

        // Create semaphore to limit queue length
        // and avoid huge memory allocation for all send and receive buffers
        unsigned int queueMaxLength = mpiCommSize;
#ifdef SEM_VALUE_MAX
        queueMaxLength = SEM_VALUE_MAX < queueMaxLength ? SEM_VALUE_MAX : queueMaxLength;
#endif
        sem_init(&masterQueue.semQueue, 0, queueMaxLength);
        pthread_create(&masterQueue.queueThread, NULL, masterQueueRun, 0);
        masterQueue.queueCreated = true;
    }
}

#ifndef MASTER_QUEUE_TEST
static void assertMPIInitialized() {
    int mpiInitialized = 0;
    MPI_Initialized(&mpiInitialized);
    assert(mpiInitialized);
}
#endif

/**
 * @brief masterQueueRun Thread entry point. This is run from initMasterQueue()
 * @param unusedArgument
 * @return 0, always
 */

static void *masterQueueRun(void *unusedArgument) {

    // dispatch queued work packages
    while(1) {
        // check if any job finished
        MPI_Status status;
        int finishedWorkerIdx = 0;

        // handle all finished jobs, if any
        while(1) {
            int dummy;
            MPI_Testany(masterQueue.numWorkers, masterQueue.recvRequests, &finishedWorkerIdx, &dummy, &status);

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
            for(int i = 0; i < masterQueue.numWorkers; ++i) {
                if(masterQueue.recvRequests[i] == MPI_REQUEST_NULL) {
                    freeWorkerIndex = i;
                    break;
                }
            }
        }

        if(freeWorkerIndex < 0) {
            // all workers are busy, wait for next one to finish
            MPI_Status status;
            MPI_Waitany(masterQueue.numWorkers, masterQueue.recvRequests, &freeWorkerIndex, &status);

            assert(freeWorkerIndex != MPI_UNDEFINED);
            receiveFinished(freeWorkerIndex, status.MPI_TAG);
        }

        masterQueueElement *currentQueueElement = getNextJob();

        if(currentQueueElement) {
            sendToWorker(freeWorkerIndex, currentQueueElement);
            masterQueue.sentJobsData[freeWorkerIndex] = currentQueueElement->data;
            free(currentQueueElement);
        }

        pthread_yield();
    };

    return 0;
}

/**
 * @brief getNextJob Pop oldest element from the queue and return.
 * @return The first queue element.
 */
static masterQueueElement *getNextJob() {

    pthread_mutex_lock(&masterQueue.mutexQueue);

    masterQueueElement *oldStart = masterQueue.queueStart;

    if(masterQueue.queueStart) {
        masterQueue.queueStart = oldStart->nextElement;
    }

    pthread_mutex_unlock(&masterQueue.mutexQueue);

    return oldStart;
}

/**
 * @brief sendToWorker Send the given work package to the given worker and track requests
 * @param workerIdx
 * @param queueElement
 */
static void sendToWorker(int workerIdx, masterQueueElement *queueElement) {
    int tag = queueElement->jobId;
    int workerRank = workerIdx + 1;
    queueData *data = queueElement->data;

#ifdef MASTER_QUEUE_H_SHOW_COMMUNICATION
    printf("\x1b[36mSent job %d to %d.\x1b[0m\n", tag, workerRank);
#endif

    MPI_Isend(data->sendBuffer, data->lenSendBuffer, MPI_BYTE, workerRank, tag, MPI_COMM_WORLD, &masterQueue.sendRequests[workerIdx]);

    MPI_Irecv(data->recvBuffer, data->lenRecvBuffer, MPI_BYTE, workerRank, tag, MPI_COMM_WORLD, &masterQueue.recvRequests[workerIdx]);

    sem_post(&masterQueue.semQueue);
}

void queueSimulation(queueData *jobData) {
    assert(masterQueue.queueCreated);

    sem_wait(&masterQueue.semQueue);

    masterQueueElement *queueElement = malloc(sizeof(*queueElement));

    queueElement->data = jobData;
    queueElement->nextElement = 0;
    masterQueueAppendElement(queueElement);
}

/**
 * @brief masterQueueAppendElement Assign job ID and append to queue.
 * @param queueElement
 */
static void masterQueueAppendElement(masterQueueElement *queueElement) {
    pthread_mutex_lock(&masterQueue.mutexQueue);

    if(masterQueue.lastJobId == INT_MAX) // Unlikely, but prevent overflow
        masterQueue.lastJobId = 0;

    queueElement->jobId = ++masterQueue.lastJobId;

    if(masterQueue.queueStart) {
        masterQueue.queueEnd->nextElement = queueElement;
    } else {
        masterQueue.queueStart = queueElement;
    }
    masterQueue.queueEnd = queueElement;

    pthread_mutex_unlock(&masterQueue.mutexQueue);
}

void terminateMasterQueue() {
    pthread_mutex_destroy(&masterQueue.mutexQueue);
    pthread_cancel(masterQueue.queueThread);
    sem_destroy(&masterQueue.semQueue);

    if(masterQueue.sentJobsData)
        free(masterQueue.sentJobsData);
    if(masterQueue.sendRequests)
        free(masterQueue.sendRequests);
    if(masterQueue.recvRequests)
        free(masterQueue.recvRequests);

    freeQueueElements();
}

/**
 * @brief freeQueueElements Deallocate queueElement memory when non-empty queue is terminated.
 * Freeing masterQueueElement.queueData memory is the users problem.
 */
static void freeQueueElements() {
    masterQueueElement *nextElement = masterQueue.queueStart;
    while(nextElement) {
        masterQueueElement *curElement = nextElement;
        nextElement = curElement->nextElement;
        free(curElement);
    }
}

/**
 * @brief receiveFinished Message received from worker, mark job as done.
 * @param workerID
 * @param jobID
 */
static void receiveFinished(int workerID, int jobID) {

#ifdef MASTER_QUEUE_H_SHOW_COMMUNICATION
    printf("\x1b[32mReceived result for job %d from %d\x1b[0m\n", jobID, workerID + 1);
#endif

    queueData *data = masterQueue.sentJobsData[workerID];
    ++(*data->jobDone);
    pthread_mutex_lock(data->jobDoneChangedMutex);
    pthread_cond_signal(data->jobDoneChangedCondition);
    pthread_mutex_unlock(data->jobDoneChangedMutex);

    masterQueue.recvRequests[workerID] = MPI_REQUEST_NULL;
}
