#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <assert.h>
#include <time.h>

#include "include/symbolic_functions.h"
#include "dataprovider.h"

#include "include/amici.h"
#include "objectivefunction.h"
#include "localoptimization.h"

void error(const char *message) { // exit?
    printf("ERROR: %s\n", message);
}

void warning(const char *message) {
    printf("WARNING: %s\n", message);
}


int main(int argc, char **argv)
{
    double initialTheta[NUM_OPTIMIZATION_PARAMS] = {0};
    getLocalOptimum(initialTheta);
}


