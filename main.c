#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <assert.h>
#include <time.h>

#include "localoptimization.h"
#include "objectivefunction.h"

int main(int argc, char **argv)
{
    // double initialTheta[NUM_OPTIMIZATION_PARAMS] = {0};
    // getLocalOptimum(initialTheta);

    UserData udata = getMyUserData();
    getLocalOptimum(udata->am_p);
}


