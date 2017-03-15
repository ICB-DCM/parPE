#include "simulationWorker.h"
#include "loadBalancerWorker.h"
#include <string.h>
#include <assert.h>
#include <include/rdata.h>
#include <include/udata.h>
#include <include/edata.h>

// TODO get rid of:
// required for datapath
#include "../objectiveFunctionBenchmarkModel/dataprovider.h"
#include "../objectiveFunctionBenchmarkModel/objectiveFunction.h"
#include "resultwriter.h"

#ifdef USE_MPE
#include <mpe.h>
extern const int mpe_event_begin_simulate, mpe_event_end_simulate;
#endif

static UserData *udata = 0;

void doWorkerWork(UserData *_udata) {
    udata = _udata;

    int workpackageLength = getLengthWorkPackageMessage(udata->am_np);
    int resultpackageLength = getLengthResultPackageMessage(udata->am_np);

    loadBalancerWorkerRun(workpackageLength, resultpackageLength, handleWorkPackage);
}


void handleWorkPackage(char *buffer, int jobId)
{
    Datapath path;
    deserializeWorkPackageMessage(buffer, udata->am_np, &path, udata->am_p, &udata->am_sensi_meth);

#if MPI_WORKER_H_VERBOSE >= 2
    printf("[%d] Received work. ", rank); printDatapath(path); fflush(stdout);
#endif

#ifdef USE_MPE
    MPE_Log_event(mpe_event_begin_simulate, jobId, "sim");
#endif

    double startTime = MPI_Wtime();

    // run simulation
    int status, iterationsUntilSteadystate = 0;
    ExpData *expData = 0; // TODO: not needed if no llh calculation here
    ReturnData *rdata = getSteadystateSolutionForExperiment(path, udata, &status, &expData, &iterationsUntilSteadystate);
    myFreeExpData(expData);

    double endTime = MPI_Wtime();
    double timeSeconds = (endTime - startTime);

#ifdef USE_MPE
    MPE_Log_event(mpe_event_end_simulate, jobId, "sim");
#endif

    char pathStrBuf[100];
    sprintDatapath(pathStrBuf, path);
    logmessage(LOGLVL_DEBUG, "Result for %s: %e  (%d) (%.2fs)", pathStrBuf, rdata->am_llhdata[0], status, timeSeconds);

    // assert(status == 0);
    // TODO write simulation results (set jobid as attribute)

    // TODO save Y
    logSimulation(path, rdata->am_llhdata[0], rdata->am_sllhdata, timeSeconds, udata->am_np, udata->am_nx, rdata->am_xdata, rdata->am_sxdata, rdata->am_ydata, jobId, iterationsUntilSteadystate);

    // Prepare result message
    resultPackageMessage result;
    result.llh = rdata->am_llhdata[0];
    result.sllh = rdata->am_sllhdata;
    result.status = status;
    serializeResultPackageMessage(result, udata->am_np, buffer);

    freeReturnData(rdata); // is referenced by workPackage, can only deallocate after serializing
}


int getLengthWorkPackageMessage(int nTheta)
{
    return sizeof(Datapath) + sizeof(double) * (nTheta) + sizeof(int);
}

void serializeWorkPackageMessage(workPackageMessage work, int nTheta, char *buffer)
{
    size_t size = 0;

    size = sizeof(work.lenData);
    memcpy(buffer, &work.lenData, size);
    buffer += size;

    size = work.lenData;
    memcpy(buffer, work.data, size);
    buffer += size;

    size = nTheta * sizeof(double);
    memcpy(buffer, work.theta, size);
    buffer += size;

    size = sizeof(int);
    memcpy(buffer, &work.sensitivityMethod, size);
    buffer += size;

}

void deserializeWorkPackageMessage(const char *msg, int nTheta, void *data, double *theta, int *sensitivityMethod)
{
    size_t size;

    int lenData = 0;
    size = sizeof(lenData);
    memcpy(&lenData, msg, size);
    msg += size;

    size = lenData;
    memcpy(data, msg, size);
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
    return sizeof(Datapath) + sizeof(double) * (nTheta + 1);
}
