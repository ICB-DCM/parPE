#include "SimulationRunner.h"
#include "simulationWorkerAmici.h"
#include <LoadBalancerMaster.h>

SimulationRunner::SimulationRunner(LoadBalancerMaster *loadBalancer)
    : loadBalancer(loadBalancer) {}

int SimulationRunner::run(
    int numJobsTotal, int lenSendBuffer,
    std::function<UserData *(int)> getUserData,
    std::function<JobIdentifier(int)> getJobIdentifier,
    std::function<int(JobData *jobs, int numJobs)> aggregate) {
    int numJobsFinished = 0;

    // TODO: allocate and free piecewise or according to max queue length
    JobData *jobs = new JobData[numJobsTotal];

    // mutex to wait for simulations to finish
    pthread_cond_t simulationsCond = PTHREAD_COND_INITIALIZER;
    pthread_mutex_t simulationsMutex = PTHREAD_MUTEX_INITIALIZER;

    for (int simulationIdx = 0; simulationIdx < numJobsTotal; ++simulationIdx) {
        // tell worker with condition to work on, for logging and reading proper
        // UserData::k
        JobIdentifier path = getJobIdentifier(simulationIdx);

        UserData *udata = getUserData(simulationIdx);

        queueSimulation(path, &jobs[simulationIdx], udata, &numJobsFinished,
                        &simulationsCond, &simulationsMutex, lenSendBuffer);
        // printf("Queued work: "); printDatapath(path);
    }

    // wait for simulations to finish
    pthread_mutex_lock(&simulationsMutex);
    while (numJobsFinished < numJobsTotal) // TODO handle finished simulations
                                           // here, don't wait for all to
                                           // complete; stop early if errors
                                           // occured
        pthread_cond_wait(&simulationsCond, &simulationsMutex);
    pthread_mutex_unlock(&simulationsMutex);
    pthread_mutex_destroy(&simulationsMutex);
    pthread_cond_destroy(&simulationsCond);

    // unpack
    int errors = aggregate(jobs, numJobsTotal);

    delete[] jobs;

    return errors;
}

void SimulationRunner::queueSimulation(JobIdentifier path, JobData *d,
                                       UserData *udata, int *jobDone,
                                       pthread_cond_t *jobDoneChangedCondition,
                                       pthread_mutex_t *jobDoneChangedMutex,
                                       int lenSendBuffer) {

    *d = JobData(lenSendBuffer, new char[lenSendBuffer], jobDone,
                 jobDoneChangedCondition, jobDoneChangedMutex);

    JobAmiciSimulation work;
    work.data = &path;
    work.lenData = sizeof(path);
    work.sensitivityMethod = udata->sensi_meth;
    work.numSimulationParameters = udata->np;
    work.simulationParameters = udata->p;
    work.serialize(d->sendBuffer);

    loadBalancer->queueJob(d);
}
