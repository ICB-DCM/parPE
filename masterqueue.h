#ifndef MASTER_QUEUE_H
#define MASTER_QUEUE_H

#include <mpi.h>
#include <dataprovider.h>

//#define MASTER_QUEUE_H_SHOW_COMMUNICATION

// callback function to provide simulation results to
typedef void (*queueSimulationFinished_cb) (void *);

// data to be sent to workers
typedef struct queueData_tag {
    int lenSendBuffer;
    char *sendBuffer;

    int lenRecvBuffer;
    char *recvBuffer;

    bool *jobDone;
    MPI_Request *recvRequest;
} queueData;

void initMasterQueue();

void queueSimulation(queueData *jobData);

void terminateMasterQueue();
#endif
