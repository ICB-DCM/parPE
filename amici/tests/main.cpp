#include <bits/stl_tree.h>

#include "CppUTest/CommandLineTestRunner.h"
#include "CppUTest/TestHarness.h"

#include <stdlib.h>
#include <time.h>

class Model;
Model *getModel() { return NULL; }

int main(int argc, char **argv) {
    srand(time(NULL));

    return CommandLineTestRunner::RunAllTests(argc, argv);
}
