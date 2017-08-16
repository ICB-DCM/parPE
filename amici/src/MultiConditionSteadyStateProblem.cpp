#include "MultiConditionSteadyStateProblem.h"

MultiConditionSteadyStateProblem::MultiConditionSteadyStateProblem() {}

ReturnData *MultiConditionSteadyStateProblem::runAndLogSimulation(
    UserData *udata, MultiConditionDataProvider *dataProvider,
    JobIdentifier path, int jobId, int *status) {
    double startTime = MPI_Wtime();

    // run simulation
    int iterationsUntilSteadystate = 0;

    ExpData *edata =
        dataProvider->getExperimentalDataForExperimentAndUpdateUserData(
            path.idxConditions, udata);

    if (edata == NULL) {
        logmessage(LOGLVL_CRITICAL,
                   "Failed to get experiment data. Check data file. Aborting.");
        abort();
    }

    ReturnData *rdata = SteadystateSimulator::getSteadystateSolution(
        udata, edata, status, &iterationsUntilSteadystate);

    freeExpData(edata);

    double endTime = MPI_Wtime();
    double timeSeconds = (endTime - startTime);

    if (isnan(rdata->llh[0]) && *status == 0)
        *status = -1;

    char pathStrBuf[100];
    sprintJobIdentifier(pathStrBuf, path);
    logmessage(LOGLVL_DEBUG, "Result for %s: %e  (%d) (%.2fs)", pathStrBuf,
               rdata->llh[0], *status, timeSeconds);

    // check for NaNs
    if (udata->sensi >= AMI_SENSI_ORDER_FIRST)
        for (int i = 0; i < udata->np; ++i)
            if (isnan(rdata->sllh[i]))
                logmessage(LOGLVL_DEBUG, "Result for %s: contains NaN at %d",
                           pathStrBuf, i);
            else if (isinf(rdata->sllh[i]))
                logmessage(LOGLVL_DEBUG, "Result for %s: contains Inf at %d",
                           pathStrBuf, i);

    // assert(*status == 0);
    // TODO write simulation results (set jobid as attribute)

    // TODO save Y
    logSimulation(path, udata->p, rdata->llh[0], rdata->sllh, timeSeconds,
                  udata->np, udata->nx, rdata->x, rdata->sx, rdata->y, jobId,
                  iterationsUntilSteadystate);

    return rdata;
}
