#include "simulationRunner.h"
#include <LoadBalancerMaster.h>

#include <omp.h>


namespace parpe {

SimulationRunnerSimple::SimulationRunnerSimple(
        std::vector<double> const& optimizationParameters,
        amici::AMICI_sensi_order sensitivityOrder,
        std::vector<int> const& conditionIndices,
        SimulationRunnerSimple::callbackJobFinishedType callbackJobFinished,
        SimulationRunnerSimple::callbackAllFinishedType aggregate)
    : optimizationParameters(optimizationParameters),
    sensitivityOrder(sensitivityOrder),
    conditionIndices(conditionIndices),
    callbackJobFinished(callbackJobFinished),
    aggregate(aggregate)
{

}


int SimulationRunnerSimple::runDistributedMemory(LoadBalancerMaster *loadBalancer, const int maxSimulationsPerPackage)
{
    int numJobsFinished = 0;

    // mutex to wait for simulations to finish
    pthread_cond_t simulationsCond = PTHREAD_COND_INITIALIZER;
    pthread_mutex_t simulationsMutex = PTHREAD_MUTEX_INITIALIZER;

    int numJobsTotal = std::ceil(static_cast<double>(conditionIndices.size()) / maxSimulationsPerPackage);
    std::vector<JobData> jobs {static_cast<unsigned int>(numJobsTotal)};

    int numConditionsSent = 0;
    for (int jobIdx = 0; jobIdx < numJobsTotal; ++jobIdx) {

        int simulationsLeft = conditionIndices.size() - numConditionsSent;
        int simulationsCurrentPackage = std::min(simulationsLeft, maxSimulationsPerPackage);
        auto currentConditions = std::vector<int>(&conditionIndices[numConditionsSent],
                         &conditionIndices[numConditionsSent + simulationsCurrentPackage]);

        queueSimulation(loadBalancer, &jobs[jobIdx],
                        &numJobsFinished, &simulationsCond, &simulationsMutex,
                        jobIdx, optimizationParameters, sensitivityOrder, currentConditions);

        numConditionsSent += simulationsCurrentPackage;
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

int SimulationRunnerSimple::runSharedMemory(LoadBalancerWorker::messageHandlerFunc messageHandler, bool sequential)
{
    std::vector<JobData> jobs {static_cast<unsigned int>(conditionIndices.size())};

    if(sequential)
        omp_set_num_threads(1);

    #pragma omp parallel for
    for (int simulationIdx = 0; simulationIdx < (signed)conditionIndices.size(); ++simulationIdx) {
        // to resuse the parallel code and for debugging we still serialze the job data here
        auto curConditionIndices = std::vector<int> {simulationIdx};
        AmiciWorkPackageSimple work {optimizationParameters, sensitivityOrder, curConditionIndices};
        auto buffer = amici::serializeToStdVec<AmiciWorkPackageSimple>(work);

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

void SimulationRunnerSimple::queueSimulation(LoadBalancerMaster *loadBalancer,
                                             JobData *d, int *jobDone,
                                             pthread_cond_t *jobDoneChangedCondition, pthread_mutex_t *jobDoneChangedMutex, int jobIdx,
                                             std::vector<double> const& optimizationParameters,
                                             amici::AMICI_sensi_order sensitivityOrder,
                                             std::vector<int> const& conditionIndices)
{
    // TODO avoid copy optimizationParameters; reuse;; for const& in work package need to split into(de)serialize
    *d = JobData(jobDone, jobDoneChangedCondition, jobDoneChangedMutex);

    AmiciWorkPackageSimple work {optimizationParameters, sensitivityOrder, conditionIndices};
    d->sendBuffer = amici::serializeToStdVec<AmiciWorkPackageSimple>(work);

    // TODO: must ignore 2nd argument for SimulationRunnerSimple
    if(callbackJobFinished)
        d->callbackJobFinished = std::bind2nd(callbackJobFinished, jobIdx);

    loadBalancer->queueJob(d);

}

void swap(SimulationRunnerSimple::AmiciResultPackageSimple &first, SimulationRunnerSimple::AmiciResultPackageSimple &second) {
    using std::swap;
    swap(first.llh, second.llh);
    swap(first.simulationTimeSeconds, second.simulationTimeSeconds);
    swap(first.gradient, second.gradient);
    swap(first.modelOutput, second.modelOutput);
    swap(first.status, second.status);
}

} // namespace parpe
