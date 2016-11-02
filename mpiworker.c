#include "mpiworker.h"
#include <string.h>
#include <mpi.h>
#include <mpe.h>
#include <assert.h>
#include "objectivefunction.h"
#include "dataprovider.h"
#include "resultwriter.h"

extern const int mpe_event_begin_simulate, mpe_event_end_simulate;

void doWorkerWork() {
    int rank, err;
    err = MPI_Comm_rank(MPI_COMM_WORLD, &rank);

#if MPI_WORKER_H_VERBOSE >= 2
    printf("[%d] Entering doWorkerWork.\n", rank); fflush(stdout);
#endif

    UserData *udata = getMyUserData();    // TODO

    int workpackageLength = getLengthWorkPackageMessage(udata->am_np);
    int resultpackageLength = getLengthResultPackageMessage(udata->am_np);
    int bufferSize = (workpackageLength > resultpackageLength) ? workpackageLength : resultpackageLength;
    char *buffer = alloca(bufferSize);

    bool terminate = false;

    while(!terminate) {

#if MPI_WORKER_H_VERBOSE >= 2
    printf("[%d] Waiting for work.\n", rank); fflush(stdout);
#endif

        // wait for receiving data to run a single simulation
        int source = 0;
        MPI_Status mpiStatus;
        err = MPI_Recv(buffer, workpackageLength, MPI_BYTE, source, MPI_ANY_TAG, MPI_COMM_WORLD, &mpiStatus);
        //printf("W%d: Received job %d\n", rank, mpiStatus.MPI_TAG);

        if(mpiStatus.MPI_TAG == MPI_TAG_EXIT_SIGNAL)
            break;

        datapath path;
        deserializeWorkPackageMessage(buffer, udata->am_np, &path, udata->am_p, &udata->am_sensi_meth);

#if MPI_WORKER_H_VERBOSE >= 2
    printf("[%d] Received work. ", rank); printDatapath(path); fflush(stdout);
#endif

        double startTime = MPI_Wtime();
        err = MPE_Log_event(mpe_event_begin_simulate, mpiStatus.MPI_TAG, "sim");
        // run simulation
        int status = 0;
        ExpData *expData = 0; // TODO: not needed if no llh calculation here
        ReturnData *rdata = getSteadystateSolutionForExperiment(path, udata, &status, &expData);
        assert(status == 0);
        // TODO write simulation results (set jobid as attribute)
        err = MPE_Log_event(mpe_event_end_simulate, mpiStatus.MPI_TAG, "sim");
        double endTime = MPI_Wtime();
        double timeSeconds = (endTime - startTime);

        logSimulation(path, rdata->am_llhdata[0], rdata->am_sllhdata, timeSeconds, udata->am_np, mpiStatus.MPI_TAG);
        // send back whatever is needed
        //reportToMaster(rdata);
        resultPackageMessage result;
        result.llh = rdata->am_llhdata[0];
        result.sllh = rdata->am_sllhdata;
        result.status = status;
        serializeResultPackageMessage(result, udata->am_np, buffer);

#if MPI_WORKER_H_VERBOSE >= 2
    printf("[%d] Simulation done, sending results (llh: %f). ", rank, result.llh); printDatapath(path); fflush(stdout);
#endif

        //printf("W%d: Sending results %d (%fs)\n", rank, mpiStatus.MPI_TAG, timeSeconds);
        MPI_Send(buffer, resultpackageLength, MPI_BYTE, 0, mpiStatus.MPI_TAG, MPI_COMM_WORLD);

        myFreeExpData(expData);
        freeReturnData(rdata);    }

    freeUserData(udata); // TODO delete free mismatch
}


// TODO Doesn't belong here
void sendTerminationSignalToAllWorkers()
{
    int commSize;
    MPI_Comm_size(MPI_COMM_WORLD, &commSize);

    MPI_Request reqs[commSize - 1];

    for(int i = 1; i < commSize; ++i) {
        reqs[i - 1] =  MPI_REQUEST_NULL;
        MPI_Isend(MPI_BOTTOM, 0, MPI_INT, i, 0, MPI_COMM_WORLD, &reqs[i - 1]);
        // printf("Sent TERM to %d\n", i);
    }

    MPI_Waitall(commSize - 1, reqs, MPI_STATUS_IGNORE);
}


void serializeWorkPackageMessage(workPackageMessage work, int nTheta, char *buffer)
{
    size_t size = 0;

    size = sizeof(datapath);
    memcpy(buffer, &work.path, size);
    buffer += size;

    size = nTheta * sizeof(double);
    memcpy(buffer, work.theta, size);
    buffer += size;

    size = sizeof(int);
    memcpy(buffer, &work.sensitivityMethod, size);
}

void deserializeWorkPackageMessage(char *msg, int nTheta, datapath *path, double *theta, int *sensitivityMethod)
{
    size_t size;

    *path = *(datapath *)msg;
    size = sizeof(datapath);
    msg += size;

    size = nTheta * sizeof(double);
    memcpy(theta, msg, size);
    msg += size;

    *sensitivityMethod = *(int *) msg;
    size = sizeof(int);
    msg += size;
}

void serializeResultPackageMessage(resultPackageMessage result, int nTheta, char *buffer)
{
    size_t size = 0;

    size = sizeof(int);
    memcpy(buffer, &result.status, size);
    buffer += size;

    size = sizeof(double);
    memcpy(buffer, &result.llh, size);
    buffer += size;

    size = nTheta * sizeof(double);
    memcpy(buffer, result.sllh, size);
    buffer += size;
}


void deserializeResultPackageMessage(char *buffer, int nTheta, int *status, double *llh, double *sllh)
{
    size_t size;

    *status = *(int *) buffer;
    size = sizeof(int);
    buffer += size;

    *llh = *(double *) buffer;
    size = sizeof(double);
    buffer += size;

    size = nTheta * sizeof(double);
    memcpy(sllh, buffer, size);
    buffer += size;
}

int getLengthWorkPackageMessage(int nTheta)
{
    return sizeof(datapath) + sizeof(double) * (nTheta) + sizeof(int);
}

int getLengthResultPackageMessage(int nTheta)
{
    return sizeof(datapath) + sizeof(double) * (nTheta + 1);
}
