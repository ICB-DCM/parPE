#include <stdlib.h>
#include "queue.h"


typedef struct QueueElement_tag {
    void *data;
    // linked list
    struct QueueElement_tag *nextElement;
} QueueElement;

typedef struct Queue_tag {
    int numElements;
    QueueElement *queueStart;
    QueueElement *queueEnd;
} Queue;

static QueueElement *queueQueueElementForData(void *data);

Queue *queueInit()
{
    Queue *queue = malloc(sizeof *queue);
    queue->numElements = 0;
    queue->queueStart = 0;
    queue->queueEnd = 0;

    return queue;
}

void queueAppend(Queue *queue, void *newElementData)
{
    QueueElement *newElement = queueQueueElementForData(newElementData);

    if(queue->queueStart) {
        queue->queueEnd->nextElement = newElement;
    } else {
        queue->queueStart = newElement;
    }

    queue->queueEnd = newElement;
    ++queue->numElements;
}

void *queuePop(Queue *queue)
{
    void *data = 0;

    if(queue->queueStart) {
        QueueElement *oldStart = queue->queueStart;
        queue->queueStart = oldStart->nextElement;
        data = oldStart->data;
        free(oldStart);
    }

    return data;
}

int queueNumberOfElements(Queue *queue)
{
    return queue->numElements;
}

static QueueElement *queueQueueElementForData(void *data) {
    QueueElement *newElement = malloc(sizeof *newElement);

    newElement->data = data;
    newElement->nextElement = 0;

    return newElement;
}

void queueDestroy(Queue *queue) {
    QueueElement *nextElement = queue->queueStart;
    while(nextElement) {
        QueueElement *curElement = nextElement;
        nextElement = curElement->nextElement;
        free(curElement);
    }
}
