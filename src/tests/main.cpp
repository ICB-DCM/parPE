#include "CppUTest/CommandLineTestRunner.h"
#include "CppUTest/TestHarness.h"

#include <time.h>
#include <stdlib.h>


TEST_GROUP(localOptimizationIpopt)
{
    void setup() {

    }

    void teardown() {

    }
};

int main(int argc, char** argv)
{
    srand(time(NULL));

    return CommandLineTestRunner::RunAllTests(argc, argv);
}
