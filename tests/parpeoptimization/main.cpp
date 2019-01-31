#include <bits/stl_tree.h>

#include "CppUTest/CommandLineTestRunner.h"
#include "CppUTest/TestHarness.h"

#include <cstdlib>
#include <ctime>

int main(int argc, char **argv) {
    srand(time(nullptr));

    // MemoryLeakWarningPlugin::turnOffNewDeleteOverloads();

    return CommandLineTestRunner::RunAllTests(argc, argv);
}
