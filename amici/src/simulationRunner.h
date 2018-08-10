#ifndef PARPE_AMICI_SIMULATIONRUNNER_H
#define PARPE_AMICI_SIMULATIONRUNNER_H

#include "multiConditionDataProvider.h" // JobIdentifier
#include <loadBalancerWorker.h>
#include <misc.h>

#include <amici/amici.h>
#include <amici/rdata.h>
#include <amici/serialization.h>

#include <functional>
#include <vector>

#include <boost/serialization/array.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/map.hpp>

#include <gsl/gsl-lite.hpp>


namespace parpe {

class JobData;
class LoadBalancerMaster;

/**
 * @brief The SimulationRunnerSimple class queues AMICI simulations, waits for the
 * results and calls a user-provided aggregation function
 */
class SimulationRunnerSimple {
  public:

    /**
     * @brief Data to be sent to a worker to run a simulation
     */
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
    int runSharedMemory(const LoadBalancerWorker::messageHandlerFunc& messageHandler, bool sequential = false);


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

    callbackJobFinishedType callbackJobFinished = nullptr;
    callbackAllFinishedType aggregate = nullptr;
    int errors = 0;
};

void swap(SimulationRunnerSimple::AmiciResultPackageSimple& first, SimulationRunnerSimple::AmiciResultPackageSimple& second);

} // namespace parpe


namespace boost {
namespace serialization {

template <class Archive>
void serialize(Archive &ar, parpe::SimulationRunnerSimple::AmiciWorkPackageSimple &u, const unsigned int version) {
    ar & u.optimizationParameters;
    ar & u.sensitivityOrder;
    ar & u.conditionIndices;
}

template <class Archive>
void serialize(Archive &ar, parpe::SimulationRunnerSimple::AmiciResultPackageSimple &u, const unsigned int version) {
    ar & u.llh;
    ar & u.simulationTimeSeconds;
    ar & u.gradient;
    ar & u.modelOutput;
    ar & u.status;   
}

} // namespace boost
} // namespace serialization

#endif // PARPE_AMICI_SIMULATIONRUNNER_H
