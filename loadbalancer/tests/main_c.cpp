#include "CppUTest/CommandLineTestRunner.h"
#include "CppUTest/TestHarness_c.h"

#include <stdlib.h>
#include <time.h>

TEST_GROUP_C_WRAPPER(queueworker){
    TEST_GROUP_C_SETUP_WRAPPER(queueworker)
        TEST_GROUP_C_TEARDOWN_WRAPPER(queueworker)};

TEST_GROUP_C_WRAPPER(queuemaster){
    TEST_GROUP_C_SETUP_WRAPPER(queuemaster)
        TEST_GROUP_C_TEARDOWN_WRAPPER(queuemaster)};
TEST_C_WRAPPER(queuemaster, test_queueinit)
TEST_C_WRAPPER(queuemaster, test_queue)
TEST_C_WRAPPER(queuemaster, test_terminateMasterQueue_noInit)
TEST_C_WRAPPER(queuemaster, test_queue_reinitialization)

TEST_GROUP_C_WRAPPER(testQueue){TEST_GROUP_C_SETUP_WRAPPER(testQueue)
                                    TEST_GROUP_C_TEARDOWN_WRAPPER(testQueue)};
TEST_C_WRAPPER(testQueue, testInitQueue)
TEST_C_WRAPPER(testQueue, testQueueAppend)
TEST_C_WRAPPER(testQueue, testQueuePop)

int main(int argc, char **argv) {
    srand(time(NULL));

    return RUN_ALL_TESTS(argc, argv);
}
