#include "CppUTest/CommandLineTestRunner.h"
#include "CppUTest/TestHarness.h"

#include <stdlib.h>
#include <time.h>

int main(int argc, char **argv) {
    srand(time(NULL));

    return CommandLineTestRunner::RunAllTests(argc, argv);
}
