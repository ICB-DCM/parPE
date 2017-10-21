#ifndef PROBLEM_H
#define PROBLEM_H

#include "LoadBalancerMaster.h"
#include "MultiConditionDataProvider.h"
#include "multiStartOptimization.h"
#include "optimizationProblem.h"
#include <LoadBalancerWorker.h>
#include <simulationWorkerAmici.h>
#include <cmath> //NAN
#include <amici.h>
#include <multiConditionProblemResultWriter.h>

namespace parpe {

class MultiConditionDataProvider;
class MultiConditionProblemResultWriter;

class MultiConditionProblem : public OptimizationProblem,
                              public LoadBalancerWorker {

  public:
    MultiConditionProblem() = default;

    MultiConditionProblem(MultiConditionDataProvider *dataProvider);

    MultiConditionProblem(MultiConditionDataProvider *dataProvider,
                          LoadBalancerMaster *loadBalancer);

    ~MultiConditionProblem() = default;

    /**
     * @brief Evaluate cost function at `optimiziationVariables`
     * @param optimiziationVariables Current parameters
     * @param objectiveFunctionValue Out: cost at `optimiziationVariables`
     * @param objectiveFunctionGradient Out: cost gradient at
     * `optimiziationVariables`
     * @return status code, non-zero on failure
     */
    virtual int
    evaluateObjectiveFunction(const double *optimiziationVariables,
                              double *objectiveFunctionValue,
                              double *objectiveFunctionGradient) override;

    virtual int evaluateObjectiveFunction(const double *optimiziationVariables,
                                          double *objectiveFunctionValue,
                                          double *objectiveFunctionGradient,
                                          int *dataIndices, int numDataIndices);

    /**
     * @brief This function is called after each iteration. See IpOpt for
     * arguments.
     * Only some are passed for CERES.
     * @param alg_mod
     * @param iter_count
     * @param obj_value
     * @param inf_pr
     * @param inf_du
     * @param mu
     * @param d_norm
     * @param regularization_size
     * @param alpha_du
     * @param alpha_pr
     * @param ls_trials
     * @return status code, non-zero to abort optimization
     */
    virtual int intermediateFunction(int alg_mod, int iter_count,
                                     double obj_value, double inf_pr,
                                     double inf_du, double mu, double d_norm,
                                     double regularization_size,
                                     double alpha_du, double alpha_pr,
                                     int ls_trials) override;

    /**
     * @brief Called after each cost function evaluation for logging results.
     * @param parameters
     * @param objectiveFunctionValue
     * @param objectiveFunctionGradient
     * @param numFunctionCalls
     * @param timeElapsed
     */
    virtual void logObjectiveFunctionEvaluation(
        const double *parameters, double objectiveFunctionValue,
        const double *objectiveFunctionGradient, int numFunctionCalls,
        double timeElapsed) override;

    /**
     * @brief Called at the end of an optimization for logging results
     * @param optimalCost
     * @param optimalParameters
     * @param masterTime
     * @param exitStatus
     */
    virtual void logOptimizerFinished(double optimalCost,
                                      const double *optimalParameters,
                                      double masterTime,
                                      int exitStatus) override;

    /**
     * @brief earlyStopping
     * @return stop the optimization run
     */
    virtual int earlyStopping();

    /**
     * @brief Is called by worker processes to run a simulation for the given
     * condition
     * @param udata UserData for simulation. Model dimensions, sensitivity
     * options and UserData::p is set, others are not.
     * E.g. UserData::k has to be update if applicable.
     * @param dataProvider
     * @param path
     * @param jobId
     * @param resultWriter
     * @return
     */
    JobResultAmiciSimulation runAndLogSimulation(amici::UserData *udata, JobIdentifier path,
                                    int jobId);

    MultiConditionDataProvider *getDataProvider();
    virtual double const*getInitialParameters(int multiStartIndex) const override;

    JobIdentifier path;

    /**
     * @brief Callback function for loadbalancer
     * @param buffer In/out: message buffer
     * @param msgSize In/out: size (bytes) of bufferobjFunVal
     * @param jobId: In: Identifier of the job (unique up to INT_MAX)
     */
    virtual void messageHandler(std::vector<char> &buffer, int jobId) override;

    virtual double getTime() const;

    std::unique_ptr<MultiConditionProblemResultWriter> resultWriter;

  protected:

    void updateUserDataCommon(const double *simulationParameters,
                              const double *objectiveFunctionGradient);

    /**
     * @brief Run AMICI simulations for conditions with the given indices
     * @param optimizationVariables
     * @param logLikelihood
     * @param objectiveFunctionGradient
     * @param dataIndices
     * @param numDataIndices
     * @return Simulation status, != 0 indicates failure
     */
    virtual int runSimulations(const double *optimizationVariables,
                               double *logLikelihood,
                               double *objectiveFunctionGradient,
                               int *dataIndices, int numDataIndices);

    /**
     * @brief Aggregates loglikelihood received from workers.
     * @param data
     * @param logLikelihood
     * @param objectiveFunctionGradient
     * @param dataIndices
     * @param numDataIndices
     * @return *Negative* log likelihood.
     */

    int aggregateLikelihood(std::vector<JobData> &data, double *logLikelihood,
                            double *objectiveFunctionGradient, int *dataIndices,
                            int numDataIndices);

    void printObjectiveFunctionFailureMessage();

    /**
     * @brief Aggregates loglikelihood gradient received from workers.
     * @param conditionIdx
     * @param simulationGradient
     * @param objectiveFunctionGradient
     * @param Gradient of the *negative* log likelihood.
     */

    void addSimulationGradientToObjectiveFunctionGradient(
        int conditionIdx, const double *simulationGradient,
        double *objectiveFunctionGradient, int numCommon);

    void
    addSimulationGradientToObjectiveFunctionGradientConditionSpecificParameters(
        const double *simulationGradient, double *objectiveFunctionGradient,
        int numCommon, int numConditionSpecificParams,
        int firstIndexOfCurrentConditionsSpecificOptimizationParameters);

    void queueSimulation(JobIdentifier path, JobData *d, int *jobDone,
                         pthread_cond_t *jobDoneChangedCondition,
                         pthread_mutex_t *jobDoneChangedMutex,
                         int lenSendBuffer);

    void setSensitivityOptions(bool sensiRequired);

    /**
     * @brief Keep information from last evaluation to avoid recomputation for
     * same parameters
     * @param optimizationParameters
     * @param objectiveFunctionValue
     * @param objectiveFunctionGradient
     */
    void
    storeCurrentFunctionEvaluation(const double *optimizationParameters,
                                   double objectiveFunctionValue,
                                   const double *objectiveFunctionGradient);

    MultiConditionDataProvider *dataProvider = nullptr;
    LoadBalancerMaster *loadBalancer = nullptr;
    amici::Model *model = nullptr;
    std::unique_ptr<amici::UserData> udata;
    amici::UserData udataOriginal; // for saving sensitivity options which are changed depending on whether gradient is needed

    // keep track of previous results to avoid re-evaluation at the same
    // parameters (using IpOpt new_x)
    std::vector<double> lastOptimizationParameters;
    std::vector<double> lastObjectiveFunctionGradient;
    double lastObjectiveFunctionValue = NAN;

    /** Time that was spent inside AMICI simulations since creating this object or last reset  */
    double simulationTimeInS = 0.0;
};

/**
 * @brief The MultiConditionProblemGeneratorForMultiStart class generates new
 * MultiConditionProblem instances with proper DataProviders for multi-start
 * optimization
 */

class MultiConditionProblemMultiStartOptimization
    : public MultiStartOptimization {
  public:
    using MultiStartOptimization::MultiStartOptimization;

    OptimizationProblem *getLocalProblemImpl(int multiStartIndex) override;

    MultiConditionDataProvider *dp = nullptr;
    OptimizationOptions options;
    MultiConditionProblemResultWriter *resultWriter = nullptr;
    amici::Model *model = nullptr;
    LoadBalancerMaster *loadBalancer = nullptr;
};

} // namespace parpe

#endif
