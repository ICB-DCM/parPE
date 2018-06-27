#ifndef SIMULATIONRUNNER_H
#define SIMULATIONRUNNER_H

#include "MultiConditionDataProvider.h" // JobIdentifier
#include <LoadBalancerWorker.h>

#include <amici/amici.h>

#include <functional>
#include <vector>

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
    int runDistributedMemory(LoadBalancerMaster *loadBalancer, const int /* TODO maxSimulationsPerPackage = 1*/);

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


class SimulationRunnerSimple {
  public:
    struct AmiciWorkPackageSimple {
        AmiciWorkPackageSimple() = default;
        std::vector<double> optimizationParameters;
        amici::AMICI_sensi_order sensitivityOrder;
        std::vector<int> conditionIndices;
        // TODO bool sendY, ...
    };

    struct AmiciResultPackageSimple {
        AmiciResultPackageSimple() = default;
        double llh;
        double simulationTimeSeconds;
        std::vector<double> gradient;
        std::vector<double> modelOutput;
        int status;
    };


//    using getUserDataType         = std::function<std::pair<std::unique_ptr<amici::Model>,std::unique_ptr<amici::Solver>> (int)>;
//    using getJobIdentifierType    = std::function<JobIdentifier(int)>;
    using callbackJobFinishedType = std::function<void(JobData*, int)>;
    using callbackAllFinishedType = std::function<int(std::vector<JobData> &)>;

    /**
     * @brief SimulationRunner
     * @param getUserData Function to provide UserData for the given simulation index. Must be provided.
     * @param getJobIdentifier Function returning a JobIdentifier object for the given simulation index. Must be provided.
     * @param callbackJobFinished Function which is called after any finished simulation.  May be nullptr.
     * @param aggregate Function which is called after all simulations are completed. May be nullptr.
     */

    SimulationRunnerSimple(const std::vector<double> &optimizationParameters,
                           amici::AMICI_sensi_order sensitivityOrder,
                           const std::vector<int> &conditionIndices,
                           //getUserDataType getUserData,
                           //getJobIdentifierType getJobIdentifier,
                           callbackJobFinishedType callbackJobFinished = nullptr,
                           callbackAllFinishedType aggregate = nullptr);

    /**
     * @brief Dispatch simulation jobs using LoadBalancerMaster
     * @param loadBalancer
     * @return
     */
    int runDistributedMemory(LoadBalancerMaster *loadBalancer, const int maxSimulationsPerPackage = 1);

    /**
     * @brief Runs simulations within the same thread. Mostly intended for
     * debugging.
     * @param messageHandler
     * @return
     */
    int runSharedMemory(LoadBalancerWorker::messageHandlerFunc messageHandler, bool sequential = false);


private:
    void queueSimulation(LoadBalancerMaster *loadBalancer,
                         JobData *d,
                         int *jobDone,
                         pthread_cond_t *jobDoneChangedCondition,
                         pthread_mutex_t *jobDoneChangedMutex,
                         int jobIdx, const std::vector<double> &optimizationParameters,
                         amici::AMICI_sensi_order sensitivityOrder, const std::vector<int> &conditionIndices);

    std::vector<double> const& optimizationParameters;
    amici::AMICI_sensi_order sensitivityOrder;
    std::vector<int> const& conditionIndices;

//    getUserDataType getUserData = nullptr;
//    getJobIdentifierType getJobIdentifier = nullptr;
    callbackJobFinishedType callbackJobFinished = nullptr;
    callbackAllFinishedType aggregate = nullptr;
    int errors = 0;


};

} // namespace parpe


namespace boost {
namespace serialization {

template <class Archive>
void serialize(Archive &ar, parpe::SimulationRunnerSimple::AmiciWorkPackageSimple &u, const unsigned int version) {
    ar &u.optimizationParameters;
    ar &u.sensitivityOrder;
    ar &u.conditionIndices;
}
template <class Archive>
void serialize(Archive &ar, parpe::SimulationRunnerSimple::AmiciResultPackageSimple &u, const unsigned int version) {
    ar &u.llh;
    ar &u.gradient;
    ar & u.status;
    ar & u.modelOutput;
}

}
}

#endif // SIMULATIONRUNNER_H
