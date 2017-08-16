#ifndef QUEUE_MASTER_H
#define QUEUE_MASTER_H

#include <mpi.h>
#include <pthread.h>

//#define MASTER_QUEUE_H_SHOW_COMMUNICATION

/** data to be sent to workers */
typedef struct JobData_tag {
    /** auto-assigned (unique number up to MAX_INT) */
    int jobId;

    /** size of data to send */
    int lenSendBuffer;
    /** data to send */
    char *sendBuffer;

    /** size of data to receive (set when job finished) */
    int lenRecvBuffer;
    /** data to receive (set when job finished) */
    char *recvBuffer;

    /** incremented by one, once the results have been received */
    int *jobDone;

    /** is signaled after jobDone has been incremented */
    pthread_cond_t *jobDoneChangedCondition;
    /** is locked to signal jobDoneChangedCondition condition  */
    pthread_mutex_t *jobDoneChangedMutex;
} JobData;

#ifdef __cplusplus
#define EXTERNC extern "C"
#else
#define EXTERNC
#endif

/**
 * @brief loadBalancerStartMaster Intialize queue. There can only be one queue.
 * Repeated calls won't do anything. MPI_Init has to be called before.
 */

EXTERNC void loadBalancerStartMaster();

/**
 * @brief loadBalancerQueueJob Append work to queue.
 * @param data
 */

EXTERNC void loadBalancerQueueJob(JobData *data);

/**
 * @brief loadBalancerTerminate Cancel the queue thread and clean up. Do not
 * wait for finish.
 */

EXTERNC void loadBalancerTerminate();

EXTERNC void sendTerminationSignalToAllWorkers();

EXTERNC JobData initJobData(int lenSendBuffer, char *sendBuffer, int *jobDone,
                            pthread_cond_t *jobDoneChangedCondition,
                            pthread_mutex_t *jobDoneChangedMutex);

#endif
