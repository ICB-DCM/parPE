#ifndef PARPE_AMICI_SIMULATIONRUNNER_H
#define PARPE_AMICI_SIMULATIONRUNNER_H

#include <parpecommon/parpeConfig.h>

#ifdef PARPE_ENABLE_MPI
#include <parpeloadbalancer/loadBalancerWorker.h>
#endif

#include <parpecommon/misc.h>

#include <amici/amici.h>
#include <amici/rdata.h>
#include <amici/serialization.h>

#include <functional>
#include <vector>

#include <boost/serialization/array.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/vector.hpp>

#include <gsl/gsl-lite.hpp>

namespace parpe {

class JobData;
#ifdef PARPE_ENABLE_MPI
class LoadBalancerMaster;
#else
// Workaround to allow building without MPI. Should be cleaned up.
using LoadBalancerMaster = int;
#endif

/**
 * @brief The AmiciSimulationRunner class queues AMICI simulations, waits for
 * the results and calls a user-provided aggregation function
 */
class AmiciSimulationRunner
{
  public:
    using messageHandlerFunc =
      std::function<void(std::vector<char>& buffer, int jobId)>;

    /**
     * @brief Data to be sent to a worker to run a simulation
     */
    struct AmiciWorkPackageSimple
    {
        AmiciWorkPackageSimple() = default;
        std::vector<double> optimizationParameters;
        amici::SensitivityOrder sensitivityOrder;
        std::vector<int> conditionIndices;
        std::string logPrefix;
        // TODO bool sendY, ...
    };

    /**
     * @brief Result from a single AMICI simulation
     */
    struct AmiciResultPackageSimple
    {
        AmiciResultPackageSimple() = default;
        double llh;
        double simulationTimeSeconds;
        std::vector<double> gradient;
        std::vector<double> modelOutput;
        std::vector<double> modelStates;
        int status;
    };

    /** Type of function be called after a single job finished  */
    using callbackJobFinishedType = std::function<void(JobData*, int)>;

    /** Type of function be called after all jobs are finished  */
    using callbackAllFinishedType = std::function<int(std::vector<JobData>&)>;

    /**
     * @brief SimulationRunner
     * @param optimizationParameters
     * @param sensitivityOrder
     * @param conditionIndices
     * @param callbackJobFinished Function which is called after any finished
     * simulation.  May be nullptr.
     * @param aggregate Function which is called after all simulations are
     * completed. May be nullptr.
     * @param logPrefix
     */
    AmiciSimulationRunner(const std::vector<double>& optimizationParameters,
                          amici::SensitivityOrder sensitivityOrder,
                          const std::vector<int>& conditionIndices,
                          callbackJobFinishedType callbackJobFinished = nullptr,
                          callbackAllFinishedType aggregate = nullptr,
                          std::string logPrefix = "");

    AmiciSimulationRunner(AmiciSimulationRunner const& other) = delete;

#ifdef PARPE_ENABLE_MPI
    /**
     * @brief Dispatch simulation jobs using LoadBalancerMaster
     * @param loadBalancer
     * @param maxSimulationsPerPackage
     * @return
     */
    int runDistributedMemory(LoadBalancerMaster* loadBalancer,
                             const int maxSimulationsPerPackage = 1);
#endif

    /**
     * @brief Runs simulations within the same thread. Mostly intended for
     * debugging.
     * @param messageHandler
     * @param sequential Run sequential (not in parallel)
     * @return
     */
    int runSharedMemory(const messageHandlerFunc& messageHandler,
                        bool sequential = false);

  private:
#ifdef PARPE_ENABLE_MPI
    void queueSimulation(LoadBalancerMaster* loadBalancer,
                         JobData* d,
                         int* jobDone,
                         pthread_cond_t* jobDoneChangedCondition,
                         pthread_mutex_t* jobDoneChangedMutex,
                         int jobIdx,
                         const std::vector<double>& optimizationParameters,
                         amici::SensitivityOrder sensitivityOrder,
                         const std::vector<int>& conditionIndices);
#endif

    std::vector<double> const& optimization_parameters_;
    amici::SensitivityOrder sensitivity_order_;
    std::vector<int> const& condition_indices_;

    callbackJobFinishedType callback_job_finished_ = nullptr;
    callbackAllFinishedType aggregate_ = nullptr;
    int errors_ = 0;
    std::string log_prefix_;
};

void
swap(AmiciSimulationRunner::AmiciResultPackageSimple& first,
     AmiciSimulationRunner::AmiciResultPackageSimple& second) noexcept;

bool
operator==(AmiciSimulationRunner::AmiciResultPackageSimple const& lhs,
           AmiciSimulationRunner::AmiciResultPackageSimple const& rhs);

} // namespace parpe

namespace boost {
namespace serialization {

template<class Archive>
void
serialize(Archive& ar,
          parpe::AmiciSimulationRunner::AmiciWorkPackageSimple& u,
          const unsigned int /*version*/)
{
    ar& u.optimizationParameters;
    ar& u.sensitivityOrder;
    ar& u.conditionIndices;
    ar& u.logPrefix;
}

template<class Archive>
void
serialize(Archive& ar,
          parpe::AmiciSimulationRunner::AmiciResultPackageSimple& u,
          const unsigned int /*version*/)
{
    ar& u.llh;
    ar& u.simulationTimeSeconds;
    ar& u.gradient;
    ar& u.modelOutput;
    ar& u.modelStates;
    ar& u.status;
}

} // namespace boost
} // namespace serialization

#endif // PARPE_AMICI_SIMULATIONRUNNER_H
