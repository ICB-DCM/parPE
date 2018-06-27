#ifndef SIMULATIONRUNNER_H
#define SIMULATIONRUNNER_H

#include "MultiConditionDataProvider.h" // JobIdentifier
#include <LoadBalancerWorker.h>

#include <amici/amici.h>
#include <amici/rdata.h>
#include <amici/serialization.h>

#include <functional>
#include <vector>

#include <misc.h>


#include <boost/serialization/array.hpp>
#include <boost/serialization/vector.hpp>

#include <gsl/gsl-lite.hpp>


namespace parpe {

class JobData;
class LoadBalancerMaster;

/**
 * @brief The JobAmiciSimulation struct contains data to be sent to a worker to run a simulation.
 */
template<typename USERDATA>
struct JobAmiciSimulation {
  public:
    JobAmiciSimulation(amici::Solver *solver, amici::Model *model, USERDATA *data)
        : solver(solver), model(model), data(data)
    {}

    std::vector<char> serialize() {
        std::string serialized;
        ::boost::iostreams::back_insert_device<std::string> inserter(serialized);
        ::boost::iostreams::stream<::boost::iostreams::back_insert_device<std::string>>
            s(inserter);
        ::boost::archive::binary_oarchive oar(s);
        // older version boost::serialization can only serialize lvalues?!
        auto parameters = model->getParameters();
        auto sensitivityMethod = solver->getSensitivityMethod();
        oar << *data
            << parameters
            << sensitivityMethod;

        s.flush();

        return std::vector<char>(serialized.begin(), serialized.end());
    }

    void deserialize(gsl::span<const char> buffer) {
        ::boost::iostreams::basic_array_source<char> device(buffer.data(), buffer.size());
        ::boost::iostreams::stream<::boost::iostreams::basic_array_source<char>> s(
            device);
        ::boost::archive::binary_iarchive iar(s);
        iar >> *data;

        std::vector<double> parameters;
        iar >> parameters;
        model->setParameters(parameters);

        int sensitivityMethod;
        iar >> sensitivityMethod;
        solver->setSensitivityMethod(static_cast<amici::AMICI_sensi_meth>(sensitivityMethod));
        solver->setSensitivityOrder(sensitivityMethod > 0 ?
                                        amici::AMICI_SENSI_ORDER_FIRST
                                      : amici::AMICI_SENSI_ORDER_NONE);
    }

    amici::Solver *solver = nullptr;
    amici::Model *model = nullptr;

    /** Simulation data or dataset Id */
    USERDATA *data = nullptr;

};


class JobResultAmiciSimulation {
public:
    JobResultAmiciSimulation() = default;
    JobResultAmiciSimulation(std::unique_ptr<amici::ReturnData> rdata, double simulationTimeInSec)
        : rdata(std::move(rdata)),
          simulationTimeInSec(simulationTimeInSec)
    {
        if(rdata)
            status = rdata->status;
    }

    ~JobResultAmiciSimulation() = default;

    friend void swap(JobResultAmiciSimulation& first, JobResultAmiciSimulation& second) {
        using std::swap;
        swap(first.status, second.status);
        swap(first.simulationTimeInSec, second.simulationTimeInSec);
        swap(first.rdata, second.rdata);
    }

    JobResultAmiciSimulation ( JobResultAmiciSimulation && other) {
        swap(*this, other);
    }

    /** simulation return status */
    int status = 1;

    /** log likelihood */
    std::unique_ptr<amici::ReturnData> rdata;

    double simulationTimeInSec = -1;
};
} // namespace parpe


namespace boost {
namespace serialization {

template <class Archive>
void serialize(Archive &ar, parpe::JobResultAmiciSimulation &d, const unsigned int version) {
    ar & d.status;
    if (Archive::is_loading::value) {
        d.rdata = std::make_unique<amici::ReturnData>();
    }
    ar & *d.rdata;
    ar & d.simulationTimeInSec;
}

} // namespace serialization
} // namespace boost


namespace parpe {


/**
 * @brief The SimulationRunner class queues AMICI simulations, waits for the
 * results and calls a user-provided aggregation function
 */
class SimulationRunner {
  public:
    using getModelAndSolver       = std::function<std::pair<std::unique_ptr<amici::Model>,std::unique_ptr<amici::Solver>> (int)>;
    using getJobIdentifierType    = std::function<JobIdentifier(int)>;
    using callbackJobFinishedType = std::function<void(JobData*, int)>;
    using callbackAllFinishedType = std::function<int(std::vector<JobData> &)>;

    /**
     * @brief SimulationRunner
     * @param numJobsTotal Total number of jobs to be run
     * @param getModelAndSolver Function to provide amici::Model and amici::Solver instances
     *  for the given simulation index. Must be provided.
     * @param getJobIdentifier Function returning a JobIdentifier object for the given simulation index. Must be provided.
     * @param callbackJobFinished Function which is called after any finished simulation.  May be nullptr.
     * @param aggregate Function which is called after all simulations are completed. May be nullptr.
     */    

    SimulationRunner(int numJobsTotal,
                     getModelAndSolver getUserData,
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
    getModelAndSolver getUserData = nullptr;
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
    ar & u.optimizationParameters;
    ar & u.sensitivityOrder;
    ar & u.conditionIndices;
}
template <class Archive>
void serialize(Archive &ar, parpe::SimulationRunnerSimple::AmiciResultPackageSimple &u, const unsigned int version) {
    ar & u.llh;
    ar & u.gradient;
    ar & u.status;
    ar & u.modelOutput;
}

}
}

#endif // SIMULATIONRUNNER_H
