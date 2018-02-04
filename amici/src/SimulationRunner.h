#ifndef SIMULATIONRUNNER_H
#define SIMULATIONRUNNER_H

#include <MultiConditionDataProvider.h> // JobIdentifier
#include <functional>
#include <vector>
#include <amici.h>
#include <LoadBalancerWorker.h>

namespace parpe {

class JobData;
class LoadBalancerMaster;

/**
 * @brief The SimulationRunner class queues AMICI simulations, waits for the
 * results and calls a user-provided aggregation function
 */
class SimulationRunner {
  public:
    using getUserDataType         = std::function<std::pair<std::unique_ptr<amici::Model>,std::unique_ptr<amici::Solver>> (int)>;
    using getJobIdentifierType    = std::function<JobIdentifier(int)>;
    using callbackJobFinishedType = std::function<void(JobData*, int)>;
    using callbackAllFinishedType = std::function<int(std::vector<JobData> &)>;

    /**
     * @brief SimulationRunner
     * @param getUserData Function to provide UserData for the given simulation index. Must be provided.
     * @param getJobIdentifier Function returning a JobIdentifier object for the given simulation index. Must be provided.
     * @param callbackJobFinished Function which is called after any finished simulation.  May be nullptr.
     * @param aggregate Function which is called after all simulations are completed. May be nullptr.
     */    

    SimulationRunner(int numJobsTotal,
                     getUserDataType getUserData,
                     getJobIdentifierType getJobIdentifier,
                     callbackJobFinishedType callbackJobFinished = nullptr,
                     callbackAllFinishedType aggregate = nullptr);

    /**
     * @brief Dispatch simulation jobs using LoadBalancerMaster
     * @param loadBalancer
     * @return
     */
    int runDistributedMemory(LoadBalancerMaster *loadBalancer);

    /**
     * @brief Runs simulations within the same thread. Mostly intended for
     * debugging.
     * @param messageHandler
     * @return
     */
    int runSharedMemory(LoadBalancerWorker::messageHandlerFunc messageHandler, bool sequential = false);


private:
    void queueSimulation(LoadBalancerMaster *loadBalancer,
                         JobIdentifier path,
                         JobData *d,
                         amici::Solver *solver, amici::Model *model,
                         int *jobDone,
                         pthread_cond_t *jobDoneChangedCondition,
                         pthread_mutex_t *jobDoneChangedMutex,
                         int simulationIdx);

    int numJobsTotal = 0;
    getUserDataType getUserData = nullptr;
    getJobIdentifierType getJobIdentifier = nullptr;
    callbackJobFinishedType callbackJobFinished = nullptr;
    callbackAllFinishedType aggregate = nullptr;
    int errors = 0;

};

} // namespace parpe

#endif // SIMULATIONRUNNER_H
