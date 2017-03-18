#ifndef QUEUE_WORKER_H
#define QUEUE_WORKER_H

#define MPI_TAG_EXIT_SIGNAL 0
#define QUEUE_WORKER_H_VERBOSE 0

/**
 * messageHandler is called by runQueueWorker when a message is received. The message is contained in buffer.
 * jobId is a message identifier, unique over the range of MAX_INT messages.
 */
typedef void (messageHandlerFp)(char* buffer, int jobId, void *userData);

void loadBalancerWorkerRun(int inMsgSize, int outMesgSize, messageHandlerFp msgHandler, void *userData);

#endif
