#ifndef SIMULATIONRUNNER_H
#define SIMULATIONRUNNER_H

#include <MultiConditionDataProvider.h> // JobIdentifier
#include <functional>
#include <vector>
#include <amici.h>

namespace parpe {

class JobData;
class LoadBalancerMaster;

/**
 * @brief The SimulationRunner class queues AMICI simulations, waits for the
 * results and calls a user-provided aggregation function
 */
class SimulationRunner {
  public:
    using getUserDataType        = std::function<amici::UserData (int)>;
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
    SimulationRunner(getUserDataType getUserData,
                     getJobIdentifierType getJobIdentifier,
                     callbackJobFinishedType callbackJobFinished,
                     callbackAllFinishedType aggregate);

    /**
     * @brief Dispatch simulation jobs using LoadBalancerMaster
     * @param numJobsTotal
     * @param lenSendBuffer
     * @param loadBalancer
     * @return
     */
    int runMPI(int numJobsTotal,
            LoadBalancerMaster *loadBalancer);

    /**
     * @brief Runs simulations within the same thread. Mostly intended for
     * debugging.
     * @param numJobsTotal
     * @param lenSendBuffer
     * @param messageHandler
     * @return
     */
    int runSerial(int numJobsTotal,
                  std::function<void(std::vector<char> &, int)> messageHandler);

    /**
     * @brief runSharedMemoryParallel
     * @param numJobsTotal
     * @param lenSendBuffer
     * @return
     */
    int runSharedMemoryParallel(int numJobsTotal, int lenSendBuffer);

    void queueSimulation(LoadBalancerMaster *loadBalancer, JobIdentifier path,
                         JobData *d, amici::UserData *udata, int *jobDone,
                         pthread_cond_t *jobDoneChangedCondition,
                         pthread_mutex_t *jobDoneChangedMutex, int simulationIdx);

  private:
    getUserDataType getUserData = nullptr;
    getJobIdentifierType getJobIdentifier = nullptr;
    callbackJobFinishedType callbackJobFinished = nullptr;
    callbackAllFinishedType aggregate = nullptr;
    int errors = 0;
};

} // namespace parpe

#endif // SIMULATIONRUNNER_H
