#include "CppUTest/TestHarness_c.h"
#include "CppUTestExt/MockSupport_c.h"

#include "queue.h"

TEST_GROUP_C_SETUP(testQueue) {
}

TEST_GROUP_C_TEARDOWN(testQueue) {
}

TEST_C(testQueue, testInitQueue) {
    Queue *queue = queueInit();

    CHECK_EQUAL_C_INT(0, queueNumberOfElements(queue));
    CHECK_EQUAL_C_POINTER(0, queuePop(queue));

    queueDestroy(queue, false);
}

TEST_C(testQueue, testQueueAppend) {
    Queue *queue = queueInit();

    queueAppend(queue, 0);
    CHECK_EQUAL_C_INT(1, queueNumberOfElements(queue));

    queueDestroy(queue, false);
}

TEST_C(testQueue, testQueuePop) {
    Queue *queue = queueInit();

    int i = 1;
    queueAppend(queue, &i);
    CHECK_EQUAL_C_INT(1, queueNumberOfElements(queue));

    int j = 1;
    queueAppend(queue, &j);
    CHECK_EQUAL_C_INT(2, queueNumberOfElements(queue));

    CHECK_EQUAL_C_POINTER(&i, queuePop(queue));
    CHECK_EQUAL_C_POINTER(&j, queuePop(queue));

    queueDestroy(queue, false);
}
