#ifndef QUEUE_H
#define QUEUE_H

#include <stdbool.h>

typedef void (*queueMapFP)(void *);

typedef struct Queue_tag Queue;

Queue *queueInit();

void queueAppend(Queue *queue, void *newElementData);

void *queuePop(Queue *queue);

int queueNumberOfElements(Queue *queue);

/**
 * @brief queueDestroy Deallocate queue memory.
 * Freeing memory of user-queued data is the users problem.
 */

void queueDestroy(Queue *queue, queueMapFP mapFunction);

#endif // QUEUE_H
