#ifndef QUEUE_WORKER_H
#define QUEUE_WORKER_H

#define MPI_TAG_EXIT_SIGNAL 0
#define QUEUE_WORKER_H_VERBOSE 0

#ifdef __cplusplus
#define EXTERNC extern "C"
#else
#define EXTERNC
#endif

/**
 * messageHandler is called by runQueueWorker when a message is received. The
 * message is contained in buffer.
 * jobId is a message identifier, unique over the range of MAX_INT messages.
 */
typedef void(messageHandlerFp)(char **buffer, int *size, int jobId,
                               void *userData);

EXTERNC void loadBalancerWorkerRun(messageHandlerFp msgHandler, void *userData);

#endif
