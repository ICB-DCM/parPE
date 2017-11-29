#include "SimulationRunner.h"
#include "simulationWorkerAmici.h"
#include <LoadBalancerMaster.h>
#include <LoadBalancerWorker.h>

namespace parpe {

SimulationRunner::SimulationRunner(int numJobsTotal,
                                   getUserDataType getUserData,
                                   getJobIdentifierType getJobIdentifier,
                                   callbackJobFinishedType callbackJobFinished)
    : SimulationRunner(numJobsTotal, getUserData, getJobIdentifier, callbackJobFinished, nullptr) {}

SimulationRunner::SimulationRunner(int numJobsTotal,
                                   getUserDataType getUserData,
                                   getJobIdentifierType getJobIdentifier,
                                   callbackJobFinishedType callbackJobFinished,
                                   callbackAllFinishedType aggregate)
    : numJobsTotal(numJobsTotal),
      getUserData(getUserData),
      getJobIdentifier(getJobIdentifier),
      callbackJobFinished(callbackJobFinished),
      aggregate(aggregate) {}

int SimulationRunner::runMPI(LoadBalancerMaster *loadBalancer) {
    int numJobsFinished = 0;

    // TODO: allocate and free piecewise or according to max queue length
    std::vector<JobData> jobs {static_cast<unsigned int>(numJobsTotal)};

    // mutex to wait for simulations to finish
    pthread_cond_t simulationsCond = PTHREAD_COND_INITIALIZER;
    pthread_mutex_t simulationsMutex = PTHREAD_MUTEX_INITIALIZER;

    for (int simulationIdx = 0; simulationIdx < numJobsTotal; ++simulationIdx) {
        // tell worker which condition to work on, for logging and reading proper
        // UserData::k
        JobIdentifier path = getJobIdentifier(simulationIdx);

        amici::UserData udata = getUserData(simulationIdx);

        queueSimulation(loadBalancer, path, &jobs[simulationIdx], &udata,
                        &numJobsFinished, &simulationsCond, &simulationsMutex,
                        simulationIdx);
        // printf("Queued work: "); printDatapath(path);
    }

    // wait for simulations to finish
    pthread_mutex_lock(&simulationsMutex);
    while (numJobsFinished < numJobsTotal) // TODO don't wait for all to
                                           // complete; stop early if errors
                                           // occured
        pthread_cond_wait(&simulationsCond, &simulationsMutex);
    pthread_mutex_unlock(&simulationsMutex);
    pthread_mutex_destroy(&simulationsMutex);
    pthread_cond_destroy(&simulationsCond);

    // unpack
    if(aggregate)
        errors += aggregate(jobs);

    return errors;
}

int SimulationRunner::runSerial(LoadBalancerWorker::messageHandlerFunc messageHandler) {

    std::vector<JobData> jobs {static_cast<unsigned int>(numJobsTotal)};

    for (int simulationIdx = 0; simulationIdx < numJobsTotal; ++simulationIdx) {
        JobIdentifier path = getJobIdentifier(simulationIdx);

        amici::UserData udata = getUserData(simulationIdx);

        JobAmiciSimulation<JobIdentifier> work(&udata, &path);
        std::vector<char> buffer = work.serialize();
        messageHandler(buffer, simulationIdx);

        jobs[simulationIdx].recvBuffer = buffer;

        if(callbackJobFinished)
            callbackJobFinished(&jobs[simulationIdx], simulationIdx);
    }

    // unpack
    if(aggregate)
        errors = aggregate(jobs);

    return errors;
}

void SimulationRunner::queueSimulation(LoadBalancerMaster *loadBalancer,
                                       JobIdentifier path, JobData *d,
                                       amici::UserData *udata, int *jobDone,
                                       pthread_cond_t *jobDoneChangedCondition,
                                       pthread_mutex_t *jobDoneChangedMutex,
                                       int simulationIdx) {

    *d = JobData(jobDone, jobDoneChangedCondition, jobDoneChangedMutex);

    JobAmiciSimulation<JobIdentifier> work(udata, &path);
    d->sendBuffer = work.serialize();

    if(callbackJobFinished)
        d->callbackJobFinished = std::bind2nd(callbackJobFinished, simulationIdx);

    loadBalancer->queueJob(d);
}

} // namespace parpe
