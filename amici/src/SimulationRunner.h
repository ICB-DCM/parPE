#ifndef SIMULATIONRUNNER_H
#define SIMULATIONRUNNER_H

#include <MultiConditionDataProvider.h> // JobIdentifier
#include <functional>
#include <vector>
class UserData;
class ExpData;
class JobData;
class LoadBalancerMaster;

/**
 * @brief The SimulationRunner class queues AMICI simulations, waits for the
 * results and calls a user-provided aggregation function
 */
class SimulationRunner {
  public:
    SimulationRunner(std::function<UserData *(int)> getUserData,
                     std::function<JobIdentifier(int)> getJobIdentifier,
                     std::function<int(std::vector<JobData> &)> aggregate);

    /**
     * @brief Dispatch simulation jobs using LoadBalancerMaster
     * @param numJobsTotal
     * @param lenSendBuffer
     * @param loadBalancer
     * @return
     */
    int run(int numJobsTotal, int lenSendBuffer,
            LoadBalancerMaster *loadBalancer);

    /**
     * @brief Runs simulations within the same thread. Mostly intended for
     * debugging.
     * @param numJobsTotal
     * @param lenSendBuffer
     * @param messageHandler
     * @return
     */
    int runSerial(int numJobsTotal, int lenSendBuffer,
                  std::function<void(char **buffer, int *msgSize, int jobId)>
                      messageHandler);

    void queueSimulation(LoadBalancerMaster *loadBalancer, JobIdentifier path,
                         JobData *d, UserData *udata, int *jobDone,
                         pthread_cond_t *jobDoneChangedCondition,
                         pthread_mutex_t *jobDoneChangedMutex,
                         int lenSendBuffer);

  private:
    std::function<UserData *(int)> getUserData = nullptr;
    std::function<JobIdentifier(int)> getJobIdentifier = nullptr;
    std::function<int(std::vector<JobData> &jobs)> aggregate = nullptr;
};

#endif // SIMULATIONRUNNER_H
