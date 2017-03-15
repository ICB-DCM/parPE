#ifndef QUEUE_H
#define QUEUE_H

#endif // QUEUE_H

typedef struct Queue_tag Queue;

Queue *queueInit();

void queueAppend(Queue *queue, void *newElementData);

void *queuePop(Queue *queue);

int queueNumberOfElements(Queue *queue);

void queueDestroy(Queue *queue);
