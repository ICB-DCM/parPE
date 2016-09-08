#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <assert.h>
#include <time.h>

#include "include/symbolic_functions.h"
#include "wrapfunctions.h"
#include "dataprovider.h"

#include "include/amici.h"
#include "objectivefunction.h"

void error(const char *message) { // exit?
    printf("ERROR: %s\n", message);
}

void warning(const char *message) {
    printf("WARNING: %s\n", message);
}

void optimizeModel() {
    bool converged = 1;
    int curIteration = 0;
    do {
        ++curIteration;
        printf("Iteration %d\n", curIteration);
        fflush(stdout);

        double timepoints [] = {};
        double theta [] = {};
        clock_t timeBegin = clock();
        double j = evaluateObjectiveFunction(timepoints, theta,
                                             NUM_FIXED_PARAMS + NUM_STATE_VARIABLES,
                                             NUM_CELL_LINES, true);
        printf("Objective function value %e\n", j);
        clock_t timeEnd = clock();
        double timeElapsed = (double) (timeEnd - timeBegin) / CLOCKS_PER_SEC;
        printf("Elapsed time %fs\n", timeElapsed);
    } while(!converged);
}

int main(int argc, char **argv)
{
    optimizeModel();
    printf("Test\n");
}


