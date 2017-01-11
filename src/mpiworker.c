#include "mpiworker.h"
#include <string.h>
#include <mpi.h>
#include <mpe.h>
#include <assert.h>
#include <alloca.h>

#include <include/rdata.h>

#include "objectivefunction.h"
#include "dataprovider.h"
#include "resultwriter.h"
#include "misc.h"

extern const int mpe_event_begin_simulate, mpe_event_end_simulate;

static resultPackageMessage handleWorkPackage(const char *buffer, UserData *udata, ReturnData **prdata, int tag);

void doWorkerWork() {
    int rank, err;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

#if MPI_WORKER_H_VERBOSE >= 4
    printf("[%d] Entering doWorkerWork.\n", rank); fflush(stdout);
#endif

    UserData *udata = getMyUserData();    // TODO

    int workpackageLength = getLengthWorkPackageMessage(udata->am_np);
    int resultpackageLength = getLengthResultPackageMessage(udata->am_np);
    int bufferSize = (workpackageLength > resultpackageLength) ? workpackageLength : resultpackageLength;
    char *buffer = alloca(bufferSize);

    bool terminate = false;

    while(!terminate) {

#if MPI_WORKER_H_VERBOSE >= 3
    printf("[%d] Waiting for work.\n", rank); fflush(stdout);
#endif

        // wait for receiving data to run a single simulation
        int source = 0;
        MPI_Status mpiStatus;
        err = MPI_Recv(buffer, workpackageLength, MPI_BYTE, source, MPI_ANY_TAG, MPI_COMM_WORLD, &mpiStatus);
        //printf("W%d: Received job %d\n", rank, mpiStatus.MPI_TAG);
        if(err != MPI_SUCCESS)
            abort();

        if(mpiStatus.MPI_TAG == MPI_TAG_EXIT_SIGNAL)
            break;

        ReturnData *rdata;

        MPE_Log_event(mpe_event_begin_simulate, mpiStatus.MPI_TAG, "sim");
        resultPackageMessage result = handleWorkPackage(buffer, udata, &rdata, mpiStatus.MPI_TAG);
        MPE_Log_event(mpe_event_end_simulate, mpiStatus.MPI_TAG, "sim");

#if MPI_WORKER_H_VERBOSE >= 2
        printf("[%d] Simulation done, sending results (llh: %f). ", rank, result.llh); printDatapath(path); fflush(stdout);
#endif
        serializeResultPackageMessage(result, udata->am_np, buffer);

        MPI_Send(buffer, resultpackageLength, MPI_BYTE, 0, mpiStatus.MPI_TAG, MPI_COMM_WORLD);

        freeReturnData(rdata); // is referenced by workPackage, can only deallocate after sending
    }

    freeUserDataC(udata);
}

resultPackageMessage handleWorkPackage(const char *buffer, UserData *udata, ReturnData **prdata, int tag)
{
    datapath path;
    deserializeWorkPackageMessage(buffer, udata->am_np, &path, udata->am_p, &udata->am_sensi_meth);

#if MPI_WORKER_H_VERBOSE >= 2
printf("[%d] Received work. ", rank); printDatapath(path); fflush(stdout);
#endif

    double startTime = MPI_Wtime();

    // run simulation
    int status, iterationsUntilSteadystate = 0;
    ExpData *expData = 0; // TODO: not needed if no llh calculation here
    *prdata = getSteadystateSolutionForExperiment(path, udata, &status, &expData, &iterationsUntilSteadystate);
    myFreeExpData(expData);

    ReturnData *rdata = *prdata;

    double endTime = MPI_Wtime();
    double timeSeconds = (endTime - startTime);

    char pathStrBuf[100];
    sprintDatapath(pathStrBuf, path);
    logmessage(LOGLVL_DEBUG, "Result for %s: %e  (%d) (%.2fs)", pathStrBuf, rdata->am_llhdata[0], status, timeSeconds);

    // assert(status == 0);
    // TODO write simulation results (set jobid as attribute)


    // TODO save Y
    logSimulation(path, rdata->am_llhdata[0], rdata->am_sllhdata, timeSeconds, udata->am_np, udata->am_nx, rdata->am_xdata, rdata->am_sxdata, rdata->am_ydata, tag, iterationsUntilSteadystate);
    // send back whatever is needed
    //reportToMaster(rdata);
    resultPackageMessage result;
    result.llh = rdata->am_llhdata[0];
    result.sllh = rdata->am_sllhdata;
    result.status = status;

    return result;
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
    }
    logmessage(LOGLVL_INFO, "Sent termination signal to workers.");
    MPI_Waitall(commSize - 1, reqs, MPI_STATUS_IGNORE);
}

int getLengthWorkPackageMessage(int nTheta)
{
    return sizeof(datapath) + sizeof(double) * (nTheta) + sizeof(int);
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
    buffer += size;

}

void deserializeWorkPackageMessage(const char *msg, int nTheta, datapath *path, double *theta, int *sensitivityMethod)
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

int getLengthResultPackageMessage(int nTheta)
{
    return sizeof(datapath) + sizeof(double) * (nTheta + 1);
}

