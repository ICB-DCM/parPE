#include "mpiworker.h"
#include <mpi.h>
#include <mpe.h>
#include "objectivefunction.h"
#include "dataprovider.h"

extern const int mpe_event_begin_simulate, mpe_event_end_simulate;

workpackage receiveData() {
    workpackage wp;

    return wp;
}

void doWorkerWork() {
    int commSize, rank, err;
    err = MPI_Comm_size(MPI_COMM_WORLD, &commSize);
    err = MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    bool terminate = false;
    UserData *udata = getMyUserData();

    while(!terminate) {
        // wait for receiving data to run a single simulation
        // TODO: for first step, just send ids and use getsteadystatesolutionforexperi
        // workpackage wp = receiveData();
        int source = 0;
        MPI_Status mpiStatus;
        err = MPI_Recv(udata->am_p, udata->am_np, MPI_DOUBLE, source, MPI_ANY_TAG, MPI_COMM_WORLD, &mpiStatus);
        //printf("W%d: Received job %d\n", rank, mpiStatus.MPI_TAG);

        if(mpiStatus.MPI_TAG == 0)
            break;

        int experimentId = mpiStatus.MPI_TAG % 1000 - 1;
        int cellLineId = mpiStatus.MPI_TAG / 1000;

        double startTime = MPI_Wtime();
        err = MPE_Log_event(mpe_event_begin_simulate, mpiStatus.MPI_TAG, "sim");
        // codify iteration, cellline, experiment
        // use custom data types
        // need: UserData, ExpData

        // write some MPI logging to analyze load balance

        // run simulation
        int status = 0;
        //ReturnData rdata = getSteadystateSolution(wp.udata, wp.edata, &status);
        ExpData *expData;
        ReturnData *rdata = getSteadystateSolutionForExperiment(cellLineId, experimentId, udata, &status, &expData);
        // write results (later)
        // writeReturnData(rdata);
        err = MPE_Log_event(mpe_event_end_simulate, mpiStatus.MPI_TAG, "sim");
        double endTime = MPI_Wtime();
        double timeSeconds = (endTime - startTime);

        // send back whatever is needed
        //reportToMaster(rdata);
        //printf("W%d: Sending results %d (%fs)\n", rank, mpiStatus.MPI_TAG, timeSeconds);
        MPI_Send(rdata->am_llhdata, 1, MPI_DOUBLE, 0, mpiStatus.MPI_TAG, MPI_COMM_WORLD);
        MPI_Send(rdata->am_sllhdata, udata->am_np, MPI_DOUBLE, 0, mpiStatus.MPI_TAG, MPI_COMM_WORLD);

        freeExpData(expData);
        freeReturnData(rdata);
        // MPI_Send(rdata->am_sllhdata, udata->am_np, MPI_DOUBLE, 0, mpiStatus.MPI_TAG, MPI_COMM_WORLD);

        // check for termination signal
        //terminate = true;
    }

    freeUserData(udata);
}

void reportToMaster(ReturnData rdata) {
    // sent llh and sllh and?
    // MPI_Send(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,  MPI_COMM_WORLD);
}

void sendTerminationSignalToAllWorkers()
{
    int commSize;
    MPI_Comm_size(MPI_COMM_WORLD, &commSize);

    for(int i = 1; i < commSize; ++i) {
        int err = MPI_Send(MPI_BOTTOM, 0, MPI_INT, i, 0, MPI_COMM_WORLD);
        // printf("Sent TERM to %d\n", i);
    }
}
