#include <stdio.h>
#include "steadystateProblemParallel.h"
#include <mpi.h>
#include <logging.h>
#include <loadBalancerMaster.h>
#include <loadBalancerWorker.h>
#include "hdf5Misc.h"
#include<unistd.h>
#include "SteadyStateMultiConditionProblem.h"
#include "wrapfunctions.h"
#include "simulationWorkerAmici.h"
#include "multiConditionProblemResultWriter.h"
/*
 * This example demonstrates the use of the loadbalancer / queue for parallel ODE simulation.
 */

void initMPI(int *argc, char ***argv);

void messageHandler(char **buffer, int *size, int jobId, void *userData);

int main(int argc, char **argv)
{
    int status = 0;

    initMPI(&argc, &argv);
    initHDF5Mutex();

    const char *filename;
    if(argc == 2) {
        filename = argv[1];
    } else {
        logmessage(LOGLVL_CRITICAL, "Must provide input file as first and only argument to %s.", argv[0]);
        return 1;
    }


    int commSize;
    MPI_Comm_size(MPI_COMM_WORLD, &commSize);

//    if(commSize == 1) {
//        // run in serial mode
//        SteadystateProblemParallel problem = SteadystateProblemParallel();
//        status = getLocalOptimum(&problem);

//    } else
    {
        int mpiRank;
        MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

        const char *outfilename = "testResultWriter_rank%3d.h5";
        char outfilefull[200];
        sprintf(outfilefull, outfilename, mpiRank);
        initResultHDFFile(outfilefull, true);

        SteadyStateMultiConditionDataProvider dataProvider =
                SteadyStateMultiConditionDataProvider(filename);
        SteadyStateMultiConditionProblem problem(&dataProvider);


        if(mpiRank == 0) {
            loadBalancerStartMaster();

            status = getLocalOptimum(&problem);

            loadBalancerTerminate();
            sendTerminationSignalToAllWorkers();
        } else {
            loadBalancerWorkerRun(messageHandler, &problem);
        }
        closeResultHDFFile();
    }

    destroyHDF5Mutex();

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

void messageHandler(char** buffer, int *msgSize, int jobId, void *userData)
{
//    int mpiRank;
//    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
//    logmessage(LOGLVL_DEBUG, "Worker #%d: Job #%d received.", mpiRank, jobId);

    SteadyStateMultiConditionProblem *problem = (SteadyStateMultiConditionProblem *) userData;
    SteadyStateMultiConditionDataProvider *dataProvider = (SteadyStateMultiConditionDataProvider *)problem->getDataProvider();

    // unpack
    UserData udata = dataProvider->getModelDims(); // TODO get from buffer // TODO make sure this is full udata, not only model dins
    dataProvider->setupUserData(&udata);

    JobIdentifier path;
    JobAmiciSimulation::toUserData(*buffer, &udata, &path);
    free(*buffer);

    // work
    int status = 0;
    ReturnData *rdata = MultiConditionProblem::runAndLogSimulation(&udata, dataProvider, path, jobId, &status);

    // pack & cleanup
    *msgSize = JobResultAmiciSimulation::getLength(udata.np);
    *buffer = (char*) malloc(*msgSize);
    JobResultAmiciSimulation::serialize(rdata, &udata, status, *buffer);

    delete rdata;

/*

    SteadystateProblemParallel *problem = (SteadystateProblemParallel*) userData;
    UserData *udata = problem->udata;

    // unpack parameters
    int conditionIdx = (int) **buffer;
    int needGradient = (int) *(*buffer + sizeof(int));
    memcpy(udata->p, *buffer + 2 * sizeof(int), sizeof(double) * udata->np);
    free(*buffer);

    // read data for current conditions
    problem->readFixedParameters(conditionIdx);
    problem->readMeasurement(conditionIdx);
    problem->requireSensitivities(needGradient);

    // run simulation
    ReturnData *rdata = getSimulationResults(udata, problem->edata);

    // pack results
    *size = sizeof(double) * (udata->nplist + 1);
    *buffer = (char*) malloc(*size);
    double *doubleBuffer = (double*) *buffer;

    doubleBuffer[0] = rdata->llh[0];
    if(needGradient)
        for(int i = 0; i < udata->nplist; ++i)
            doubleBuffer[1 + i] = rdata->sllh[i];

    delete rdata;
    */
}
