#include <stdio.h>
#include "steadystateProblemParallel.h"
#include <mpi.h>
#include <logging.h>
#include <loadBalancerMaster.h>
#include <loadBalancerWorker.h>

#include<unistd.h>
/*
 * This example demonstrates the use of the loadbalancer / queue for parallel ODE simulation.
 */

void initMPI(int *argc, char ***argv);

void messageHandler(char **buffer, int *size, int jobId, void *userData);

int main(int argc, char **argv)
{
    int status = 0;

    initMPI(&argc, &argv);

    SteadystateProblemParallel problem = SteadystateProblemParallel(1);

    int commSize;
    MPI_Comm_size(MPI_COMM_WORLD, &commSize);

    if(commSize == 1) {
        // run in serial mode
        status = getLocalOptimum(&problem);
    } else {
        int mpiRank;
        MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
        if(mpiRank == 0) {
            loadBalancerStartMaster();

            status = getLocalOptimum(&problem);

            loadBalancerTerminate();
            sendTerminationSignalToAllWorkers();
        } else {
            loadBalancerWorkerRun(messageHandler, &problem);
        }
    }

    MPI_Finalize();

    return status;
}

void initMPI(int *argc, char ***argv) {
    int mpiErr = MPI_Init(argc, argv);
    if(mpiErr != MPI_SUCCESS) {
        logmessage(LOGLVL_CRITICAL, "Problem initializing MPI. Exiting.");
        exit(1);
    }

    int mpiRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    if(mpiRank == 0) {
        int commSize;
        MPI_Comm_size(MPI_COMM_WORLD, &commSize);

        logmessage(LOGLVL_INFO, "Running with %d MPI processes.", commSize);
    }
}

void messageHandler(char** buffer, int *size, int jobId, void *userData)
{
//    logmessage(LOGLVL_DEBUG, "Job #%d received.", jobId);

    SteadystateProblemParallel *problem = (SteadystateProblemParallel*) userData;
    UserData *udata = problem->udata;

    // unpack data
    double *doubleBuffer = (double *) *buffer;
    for(int i = 0; i < udata->am_nk; ++i)
        udata->am_k[i] = doubleBuffer[i];
    for(int i = 0; i < udata->am_np; ++i)
        udata->am_p[i] = doubleBuffer[i + udata->am_nk];
    free(*buffer);

    // run simulation
    int tmpStatus;
    ReturnData *rdata = getSimulationResults(udata, problem->edata, &tmpStatus);

    // pack results
    *size = sizeof(double) * (udata->am_nplist + 1);
    *buffer = (char*) malloc(*size);
    doubleBuffer = (double*) *buffer;

    doubleBuffer[0] = rdata->am_llhdata[0];
    for(int i = 0; i < udata->am_nplist; ++i)
        doubleBuffer[1 + i] = rdata->am_sllhdata[i];

    freeReturnData(rdata);
}
