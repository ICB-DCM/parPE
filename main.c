#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <assert.h>
#include <time.h>
#include <mpi.h>

#include "localoptimization.h"
#include "objectivefunction.h"

int main(int argc, char **argv)
{
    int commSize, rank;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &commSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    printf(" I am %d of %d.\n", rank, commSize);

    // double initialTheta[NUM_OPTIMIZATION_PARAMS] = {0};
    // getLocalOptimum(initialTheta);

    UserData udata = getMyUserData();
    getLocalOptimum(udata->am_p);

    MPI_Finalize();

}


