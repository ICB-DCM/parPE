#include "CppUTest/CommandLineTestRunner.h"
#include "CppUTest/TestHarness_c.h"

TEST_GROUP_C_WRAPPER(objectivefunction)
{
    TEST_GROUP_C_SETUP_WRAPPER(objectivefunction)
    TEST_GROUP_C_TEARDOWN_WRAPPER(objectivefunction)
};

TEST_C_WRAPPER(objectivefunction, test_reachedSteadyState)

int main(int argc, char** argv)
{
    return RUN_ALL_TESTS(argc, argv);
}
