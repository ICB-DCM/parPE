#ifndef MASTER_QUEUE_H
#define MASTER_QUEUE_H

#include <mpi.h>
#include <pthread.h>

//#define MASTER_QUEUE_H_SHOW_COMMUNICATION

/** data to be sent to workers */
typedef struct queueData_tag {
    int lenSendBuffer;
    char *sendBuffer;

    int lenRecvBuffer;
    char *recvBuffer;

    /** incremented by one, once the results have been received */
    int *jobDone;

    /** is signaled after jobDone has been incremented */
    pthread_cond_t *jobDoneChangedCondition;
    pthread_mutex_t *jobDoneChangedMutex;

    MPI_Request *recvRequest;
} queueData;

void initMasterQueue();

void queueSimulation(queueData *jobData);

void terminateMasterQueue();

#endif
