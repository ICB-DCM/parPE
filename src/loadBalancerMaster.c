#include "loadBalancerMaster.h"
#include "queue.h"

#include <pthread.h>
#include <semaphore.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <limits.h>

typedef struct LoadBalancer_tag {
    int numWorkers;
    Queue *queue;
    int lastJobId;
    // one for each worker, index is off by one from MPI rank
    // because no job is sent to master (rank 0)
    bool *workerIsBusy;
    MPI_Request *sendRequests; // TODO: option: free(sendbuffer)
    JobData **sentJobsData;
    pthread_mutex_t mutexQueue;
    sem_t semQueue;
    pthread_t queueThread;
} LoadBalancer;

#define QUEUEMASTER_QUEUE_INITIALIZER {0, NULL, 0, NULL, NULL, NULL, PTHREAD_MUTEX_INITIALIZER, {}, 0}

static LoadBalancer loadBalancer = QUEUEMASTER_QUEUE_INITIALIZER;

static void *loadBalancerRun(void *unusedArgument);

static void assertMPIInitialized();

static JobData *getNextJob();

static void sendToWorker(int workerIdx, JobData *queueElement);

static void masterQueueAppendElement(JobData *queueElement);

static int handleReply(MPI_Status *mpiStatus);

/**
 * @brief initLoadBalancerMaster Intialize load balancer.
 */
void loadBalancerStartMaster() {
    // There can only be one queue
    if(!loadBalancer.queue) {
        assertMPIInitialized();

        loadBalancer.queue = queueInit();

        int mpiCommSize;
        MPI_Comm_size(MPI_COMM_WORLD, &mpiCommSize);
        assert(mpiCommSize > 1); // crashes otherwise

        loadBalancer.numWorkers = mpiCommSize - 1;
        loadBalancer.sentJobsData = malloc(loadBalancer.numWorkers * sizeof(JobData *));
        loadBalancer.workerIsBusy = malloc(loadBalancer.numWorkers * sizeof(bool));
        memset(loadBalancer.workerIsBusy, 0, loadBalancer.numWorkers * sizeof(bool));
        loadBalancer.sendRequests = malloc(loadBalancer.numWorkers * sizeof(MPI_Request));
        for(int i = 0; i < loadBalancer.numWorkers; ++i) // have to initialize before can wait!
            loadBalancer.sendRequests[i] = MPI_REQUEST_NULL;

        // Create semaphore to limit queue length
        // and avoid huge memory allocation for all send and receive buffers
        unsigned int queueMaxLength = mpiCommSize;
#ifdef SEM_VALUE_MAX
        queueMaxLength = SEM_VALUE_MAX < queueMaxLength ? SEM_VALUE_MAX : queueMaxLength;
#endif
        sem_init(&loadBalancer.semQueue, 0, queueMaxLength);
        pthread_create(&loadBalancer.queueThread, NULL, loadBalancerRun, 0);
    }
}

#ifndef QUEUE_MASTER_TEST
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

static void *loadBalancerRun(void *unusedArgument) {

    // dispatch queued work packages
    while(1) {

        // check if any job finished
        MPI_Status status;
        int finishedWorkerIdx = -1;

        // handle all finished jobs, if any
        while(1) {
            // add cancellation point to avoid invalid reads in loadBalancer.recvRequests
            pthread_testcancel();

            int flag;
            MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &status);

            if(flag) {
                // some job is finished
                finishedWorkerIdx = handleReply(&status);
            } else {
                // there was nothing to be finished
                break;
            }
        }

        // getNextFreeWorker
        int freeWorkerIndex = finishedWorkerIdx;

        if(freeWorkerIndex < 0) {
            // no job finished recently, check free slots
            for(int i = 0; i < loadBalancer.numWorkers; ++i) {
                if(loadBalancer.workerIsBusy[i] == false) {
                    freeWorkerIndex = i;
                    break;
                }
            }
        }

        if(freeWorkerIndex < 0) {
            // add cancellation point to avoid invalid reads in loadBalancer.recvRequests
            pthread_testcancel();

            // all workers are busy, wait for next one to finish
            MPI_Status status;
            MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            freeWorkerIndex = handleReply(&status);
        }

        // found free worker, check for jobs to do
        JobData *currentQueueElement = getNextJob();

        if(currentQueueElement) {
            sendToWorker(freeWorkerIndex, currentQueueElement);
            loadBalancer.sentJobsData[freeWorkerIndex] = currentQueueElement;
        }

        pthread_yield();
    };

    return 0;
}

/**
 * @brief getNextJob Pop oldest element from the queue and return.
 * @return The first queue element.
 */
static JobData *getNextJob() {

    pthread_mutex_lock(&loadBalancer.mutexQueue);

    JobData *nextJob = (JobData *) queuePop(loadBalancer.queue);

    pthread_mutex_unlock(&loadBalancer.mutexQueue);

    return nextJob;
}

/**
 * @brief sendToWorker Send the given work package to the given worker and track requests
 * @param workerIdx
 * @param queueElement
 */
static void sendToWorker(int workerIdx, JobData *data) {
    assert(workerIdx >= 0);
    assert(workerIdx < loadBalancer.numWorkers);

    loadBalancer.workerIsBusy[workerIdx] = true;

    int tag = data->jobId;
    int workerRank = workerIdx + 1;

#ifdef MASTER_QUEUE_H_SHOW_COMMUNICATION
    printf("\x1b[36mSent job #%d to rank %d.\x1b[0m\n", tag, workerRank);
#endif

    MPI_Isend(data->sendBuffer, data->lenSendBuffer, MPI_BYTE, workerRank, tag, MPI_COMM_WORLD, &loadBalancer.sendRequests[workerIdx]);

    sem_post(&loadBalancer.semQueue);
}

/**
 * @brief masterQueueAppendElement Assign job ID and append to queue.
 * @param queueElement
 */

void loadBalancerQueueJob(JobData *data) {
    assert(loadBalancer.queue);

    sem_wait(&loadBalancer.semQueue);

    pthread_mutex_lock(&loadBalancer.mutexQueue);

    if(loadBalancer.lastJobId == INT_MAX) // Unlikely, but prevent overflow
        loadBalancer.lastJobId = 0;

    data->jobId = ++loadBalancer.lastJobId;

    queueAppend(loadBalancer.queue, data);

    pthread_mutex_unlock(&loadBalancer.mutexQueue);

}

void loadBalancerTerminate() {
    pthread_cancel(loadBalancer.queueThread);
    pthread_join(loadBalancer.queueThread, NULL);

    pthread_mutex_destroy(&loadBalancer.mutexQueue);
    sem_destroy(&loadBalancer.semQueue);

    if(loadBalancer.sentJobsData)
        free(loadBalancer.sentJobsData);
    if(loadBalancer.sendRequests)
        free(loadBalancer.sendRequests);
    if(loadBalancer.workerIsBusy)
        free(loadBalancer.workerIsBusy);

    queueDestroy(loadBalancer.queue, 0);
    loadBalancer.queue = 0;
}


/**
 * @brief receiveFinished Message received from worker, mark job as done.
 * @param workerID
 * @param jobID
 */
static int handleReply(MPI_Status *mpiStatus) {

    int workerIdx = mpiStatus->MPI_SOURCE - 1;
    JobData *data = loadBalancer.sentJobsData[workerIdx];

    // allocate memory for result
    int size;
    MPI_Get_count(mpiStatus, MPI_BYTE, &size);
    data->recvBuffer = malloc(size);

#ifdef MASTER_QUEUE_H_SHOW_COMMUNICATION
    printf("\x1b[32mReceiving result for job %d from %d (%dB)\x1b[0m\n", mpiStatus->MPI_TAG, mpiStatus->MPI_SOURCE, size);
#endif

    // receive
    MPI_Recv(data->recvBuffer, size, MPI_BYTE, mpiStatus->MPI_SOURCE, mpiStatus->MPI_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    ++(*data->jobDone);
    pthread_mutex_lock(data->jobDoneChangedMutex);
    pthread_cond_signal(data->jobDoneChangedCondition);
    pthread_mutex_unlock(data->jobDoneChangedMutex);

    // free send buffer TODO: this should be done earlier using TestAny on sendRequests
    loadBalancer.sendRequests[workerIdx] = MPI_REQUEST_NULL;
    free(data->sendBuffer);
    loadBalancer.workerIsBusy[workerIdx] = false;

#ifdef MASTER_QUEUE_H_SHOW_COMMUNICATION
    printf("\x1b[32mReceived result for job %d from %d\x1b[0m\n", mpiStatus->MPI_TAG, mpiStatus->MPI_SOURCE);
#endif

    return workerIdx;
}

void sendTerminationSignalToAllWorkers()
{
    int commSize;
    MPI_Comm_size(MPI_COMM_WORLD, &commSize);

    MPI_Request reqs[commSize - 1];

    for(int i = 1; i < commSize; ++i) {
        reqs[i - 1] =  MPI_REQUEST_NULL;
        MPI_Isend(MPI_BOTTOM, 0, MPI_INT, i, 0, MPI_COMM_WORLD, &reqs[i - 1]);
    }
    MPI_Waitall(commSize - 1, reqs, MPI_STATUS_IGNORE);
}
