#include <time.h>
#include <stdlib.h>

#include "CppUTest/CommandLineTestRunner.h"
#include "CppUTest/TestHarness_c.h"

TEST_GROUP_C_WRAPPER(objectivefunction)
{
    TEST_GROUP_C_SETUP_WRAPPER(objectivefunction)
    TEST_GROUP_C_TEARDOWN_WRAPPER(objectivefunction)
};
TEST_C_WRAPPER(objectivefunction, test_reachedSteadyState)

TEST_GROUP_C_WRAPPER(mpiworker)
{
    TEST_GROUP_C_SETUP_WRAPPER(mpiworker)
    TEST_GROUP_C_TEARDOWN_WRAPPER(mpiworker)
};
TEST_C_WRAPPER(mpiworker, test_serializeResultPackageMessage)
TEST_C_WRAPPER(mpiworker, test_serializeWorkPackageMessage)

TEST_GROUP_C_WRAPPER(masterqueue)
{
    TEST_GROUP_C_SETUP_WRAPPER(masterqueue)
    TEST_GROUP_C_TEARDOWN_WRAPPER(masterqueue)
};
TEST_C_WRAPPER(masterqueue, test_queueinit)
TEST_C_WRAPPER(masterqueue, test_queue)
TEST_C_WRAPPER(masterqueue, test_terminateMasterQueue_noInit)
TEST_C_WRAPPER(masterqueue, test_queue_reinitialization)


int main(int argc, char** argv)
{
    srand(time(NULL));

    return RUN_ALL_TESTS(argc, argv);
}
