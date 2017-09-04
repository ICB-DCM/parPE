#ifndef SIMULATIONRUNNER_H
#define SIMULATIONRUNNER_H

#include <MultiConditionDataProvider.h> // JobIdentifier
#include <functional>
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
    SimulationRunner(LoadBalancerMaster *loadBalancer);

    int run(int numJobsTotal, int lenSendBuffer,
            std::function<UserData *(int)> getUserData,
            std::function<JobIdentifier(int)> getJobIdentifier,
            std::function<int(JobData *jobs, int numJobs)> aggregate);

    void queueSimulation(JobIdentifier path, JobData *d, UserData *udata,
                         int *jobDone, pthread_cond_t *jobDoneChangedCondition,
                         pthread_mutex_t *jobDoneChangedMutex,
                         int lenSendBuffer);

  private:
    LoadBalancerMaster *loadBalancer;
};

#endif // SIMULATIONRUNNER_H
