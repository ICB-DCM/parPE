#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <assert.h>
#include <time.h>

#include <mpi.h>
#include <mpe.h>
#include "mpiworker.h"

#include "localoptimization.h"
#include "objectivefunction.h"
#include "logger.h"

int mpe_event_begin_simulate, mpe_event_end_simulate;
int mpe_event_begin_getrefs, mpe_event_end_getrefs;
int mpe_event_begin_aggregate, mpe_event_end_aggregate;

int main(int argc, char **argv)
{
    int commSize, rank, err;

    err = MPI_Init(&argc, &argv);
    err = MPI_Comm_size(MPI_COMM_WORLD, &commSize);
    err = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    assert(commSize > 1);

    err = MPE_Init_log();

    MPE_Log_get_state_eventIDs(&mpe_event_begin_simulate, &mpe_event_end_simulate);
    MPE_Log_get_state_eventIDs(&mpe_event_begin_aggregate, &mpe_event_end_aggregate);
    MPE_Log_get_state_eventIDs(&mpe_event_begin_getrefs, &mpe_event_end_getrefs);

    if(rank == 0) {
        err = MPE_Describe_state(mpe_event_begin_simulate, mpe_event_end_simulate, "simulate", "blue:gray");
        err = MPE_Describe_state(mpe_event_begin_aggregate, mpe_event_end_aggregate, "aggregate", "red:gray");
        err = MPE_Describe_state(mpe_event_begin_getrefs, mpe_event_end_getrefs, "getrefs", "green:gray");

        // double initialTheta[NUM_OPTIMIZATION_PARAMS] = {0};
        // getLocalOptimum(initialTheta);

        loggerdata datalogger = initResultHDFFile("results.h5", "/results/multistarts/0");
        UserData udata = getMyUserData();
        getLocalOptimum(udata->am_p, &datalogger);
        freeUserData(udata);
        closeResultHDFFile(datalogger);
        sendTerminationSignalToAllWorkers();
    } else {
        doWorkerWork();
    }

    err = MPE_Finish_log("mpe.log");
    err = MPI_Finalize();
}
