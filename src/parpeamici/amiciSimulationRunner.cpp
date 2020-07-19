#include <parpeamici/amiciSimulationRunner.h>

#include <parpeloadbalancer/loadBalancerMaster.h>

#if defined(_OPENMP)
#include <omp.h>
#endif

#include <utility>

// #define PARPE_SIMULATION_RUNNER_DEBUG

namespace parpe {

AmiciSimulationRunner::AmiciSimulationRunner(
  std::vector<double> const& optimizationParameters,
  amici::SensitivityOrder sensitivityOrder,
  std::vector<int> const& conditionIndices,
  AmiciSimulationRunner::callbackJobFinishedType callbackJobFinished,
  AmiciSimulationRunner::callbackAllFinishedType aggregate,
  std::string logPrefix)
  : optimization_parameters_(optimizationParameters)
  , sensitivity_order_(sensitivityOrder)
  , condition_indices_(conditionIndices)
  , callback_job_finished_(std::move(std::move(callbackJobFinished)))
  , aggregate_(std::move(std::move(aggregate)))
  , log_prefix_(std::move(logPrefix))
{}

#ifdef PARPE_ENABLE_MPI
int
AmiciSimulationRunner::runDistributedMemory(LoadBalancerMaster* loadBalancer,
                                            const int maxSimulationsPerPackage)
{
#ifdef PARPE_SIMULATION_RUNNER_DEBUG
    printf("runDistributedMemory\n");
#endif

    // mutex and condition to wait for simulations to finish
    pthread_cond_t simulationsCond = PTHREAD_COND_INITIALIZER;
    pthread_mutex_t simulationsMutex = PTHREAD_MUTEX_INITIALIZER;

    // multiple simulations may be grouped into one work package
    auto numJobsTotal = static_cast<int>(
      std::ceil(static_cast<double>(condition_indices_.size()) /
                maxSimulationsPerPackage));
    std::vector<JobData> jobs{ static_cast<decltype(jobs)::size_type>(
      numJobsTotal) };
    int numJobsFinished = 0;
    int numConditionsSent = 0;

    // prepare and queue work package
    for (int jobIdx = 0; jobIdx < numJobsTotal; ++jobIdx) {
        int simulationsLeft =
          static_cast<int>(condition_indices_.size()) - numConditionsSent;
        int simulationsCurrentPackage =
          std::min(simulationsLeft, maxSimulationsPerPackage);

        auto currentConditions = std::vector<int>(
          &condition_indices_[static_cast<std::vector<int>::size_type>(
            numConditionsSent)],
          &condition_indices_[numConditionsSent + simulationsCurrentPackage]);
        queueSimulation(loadBalancer,
                        &jobs[jobIdx],
                        &numJobsFinished,
                        &simulationsCond,
                        &simulationsMutex,
                        jobIdx,
                        optimization_parameters_,
                        sensitivity_order_,
                        currentConditions);

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
    if (aggregate_)
        errors_ += aggregate_(jobs);

    return errors_;
}
#endif

int
AmiciSimulationRunner::runSharedMemory(const messageHandlerFunc& messageHandler,
                                       bool sequential)
{
#ifdef PARPE_SIMULATION_RUNNER_DEBUG
    printf("runSharedMemory\n");
#endif

    std::vector<JobData> jobs{ static_cast<unsigned int>(
      condition_indices_.size()) };

#if defined(_OPENMP)
    if (sequential)
        omp_set_num_threads(1);

#pragma omp parallel for
#endif
    for (int simulationIdx = 0;
         simulationIdx < (signed)condition_indices_.size();
         ++simulationIdx) {
        // to reuse the parallel code and for debugging we still serialze the
        // job data here
        auto curConditionIndices = std::vector<int>{ simulationIdx };
        AmiciWorkPackageSimple work{ optimization_parameters_,
                                     sensitivity_order_,
                                     curConditionIndices,
                                     log_prefix_ };
        auto buffer = amici::serializeToStdVec<AmiciWorkPackageSimple>(work);

        messageHandler(buffer, simulationIdx);
        jobs[simulationIdx].recvBuffer = buffer;

        if (callback_job_finished_)
            callback_job_finished_(&jobs[simulationIdx], simulationIdx);
    }

    // unpack
    if (aggregate_)
        errors_ = aggregate_(jobs);

    return errors_;
}

#ifdef PARPE_ENABLE_MPI
void
AmiciSimulationRunner::queueSimulation(
  LoadBalancerMaster* loadBalancer,
  JobData* d,
  int* jobDone,
  pthread_cond_t* jobDoneChangedCondition,
  pthread_mutex_t* jobDoneChangedMutex,
  int jobIdx,
  std::vector<double> const& optimizationParameters,
  amici::SensitivityOrder sensitivityOrder,
  std::vector<int> const& conditionIndices)
{
    // TODO avoid copy optimizationParameters; reuse;; for const& in work
    // package need to split into(de)serialize
    *d = JobData(jobDone, jobDoneChangedCondition, jobDoneChangedMutex);

    AmiciWorkPackageSimple work{
        optimizationParameters, sensitivityOrder, conditionIndices, log_prefix_
    };
    d->sendBuffer = amici::serializeToStdVec<AmiciWorkPackageSimple>(work);

    // TODO: must ignore 2nd argument for SimulationRunnerSimple
    if (callback_job_finished_)
        d->callbackJobFinished = std::bind2nd(callback_job_finished_, jobIdx);

    loadBalancer->queueJob(d);
}
#endif

void
swap(AmiciSimulationRunner::AmiciResultPackageSimple& first,
     AmiciSimulationRunner::AmiciResultPackageSimple& second) noexcept
{
    using std::swap;
    swap(first.llh, second.llh);
    swap(first.simulationTimeSeconds, second.simulationTimeSeconds);
    swap(first.gradient, second.gradient);
    swap(first.modelOutput, second.modelOutput);
    swap(first.modelStates, second.modelStates);
    swap(first.modelSigmas, second.modelSigmas);
    swap(first.status, second.status);
}

bool
operator==(const AmiciSimulationRunner::AmiciResultPackageSimple& lhs,
           const AmiciSimulationRunner::AmiciResultPackageSimple& rhs)
{
    return lhs.llh == rhs.llh && lhs.status == rhs.status &&
           lhs.gradient == rhs.gradient && lhs.modelOutput == rhs.modelOutput &&
           lhs.modelStates == rhs.modelStates &&
           lhs.modelSigmas == rhs.modelSigmas &&
           lhs.simulationTimeSeconds == rhs.simulationTimeSeconds;
}

} // namespace parpe
