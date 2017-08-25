#include "hdf5Misc.h"
#include "steadystateProblemParallel.h"
#include <amici_model.h>
#include <LoadBalancerMaster.h>
#include <loadBalancerWorker.h>
#include <logging.h>
#include <mpi.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
/*
 * This example demonstrates the use of the loadbalancer / queue for parallel
 * ODE simulation.
 */

void initMPI(int *argc, char ***argv);

void messageHandler(char **buffer, int *size, int jobId, void *userData);

int main(int argc, char **argv) {
    int status = 0;

    initMPI(&argc, &argv);
    initHDF5Mutex();

    int commSize;
    MPI_Comm_size(MPI_COMM_WORLD, &commSize);

    if (commSize == 1) {
        // run in serial mode
        SteadystateProblemParallel problem = SteadystateProblemParallel(NULL);
        status = getLocalOptimum(&problem);

    } else {
        SteadystateProblemParallel problem = SteadystateProblemParallel(NULL);

        int mpiRank;
        MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

        if (mpiRank == 0) {
            LoadBalancerMaster lbm;
            problem.loadBalancer = &lbm;

            lbm.run();

            status = getLocalOptimum(&problem);

            lbm.terminate();
            lbm.sendTerminationSignalToAllWorkers();
        } else {
            loadBalancerWorkerRun(messageHandler, &problem);
        }
    }

    destroyHDF5Mutex();

    MPI_Finalize();

    return status;
}

void initMPI(int *argc, char ***argv) {
    int mpiErr = MPI_Init(argc, argv);
    if (mpiErr != MPI_SUCCESS) {
        logmessage(LOGLVL_CRITICAL, "Problem initializing MPI. Exiting.");
        exit(1);
    }

    int mpiRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    if (mpiRank == 0) {
        int commSize;
        MPI_Comm_size(MPI_COMM_WORLD, &commSize);

        logmessage(LOGLVL_INFO, "Running with %d MPI processes.", commSize);
    }
}

void messageHandler(char **buffer, int *size, int jobId, void *userData) {
    //    int mpiRank;
    //    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    //    logmessage(LOGLVL_DEBUG, "Worker #%d: Job #%d received.", mpiRank,
    //    jobId);

    SteadystateProblemParallel *problem =
        (SteadystateProblemParallel *)userData;
    UserData *udata = problem->udata;
    Model *model = problem->model;

    // unpack parameters
    int conditionIdx = (int)**buffer;
    int needGradient = (int)*(*buffer + sizeof(int));
    memcpy(udata->p, *buffer + 2 * sizeof(int), sizeof(double) * model->np);
    free(*buffer);

    // read data for current conditions
    problem->readFixedParameters(conditionIdx);
    problem->readMeasurement(conditionIdx);
    problem->requireSensitivities(needGradient);

    // run simulation
    ReturnData *rdata = getSimulationResults(model, udata, problem->edata);
    // printf("Result for %d: %f\n", conditionIdx, *rdata->llh);
    // pack results
    *size = sizeof(double) * (udata->nplist + 1);
    *buffer = (char *)malloc(*size);
    double *doubleBuffer = (double *)*buffer;

    doubleBuffer[0] = rdata->llh[0];
    if (needGradient)
        for (int i = 0; i < udata->nplist; ++i)
            doubleBuffer[1 + i] = rdata->sllh[i];

    delete rdata;
}
