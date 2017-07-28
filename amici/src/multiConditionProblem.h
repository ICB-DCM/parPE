#ifndef PROBLEM_H
#define PROBLEM_H

#include "optimizationProblem.h"
#include "multiStartOptimization.h"
#include "multiConditionProblemResultWriter.h"
#include "MultiConditionDataProvider.h"
#include "amici_interface_cpp.h"
#include "loadBalancerMaster.h"
#include "rdata.h"

class MultiConditionProblem : public OptimizationProblem {

public:
    MultiConditionProblem();

    MultiConditionProblem(MultiConditionDataProvider *dataProvider);

    /**
     * @brief Evaluate cost function at `optimiziationVariables`
     * @param optimiziationVariables Current parameters
     * @param objectiveFunctionValue Out: cost at `optimiziationVariables`
     * @param objectiveFunctionGradient Out: cost gradient at `optimiziationVariables`
     * @return status code, non-zero on failure
     */
    virtual int evaluateObjectiveFunction(const double *optimiziationVariables, double *objectiveFunctionValue, double *objectiveFunctionGradient);

    virtual int evaluateObjectiveFunction(const double *optimiziationVariables, double *objectiveFunctionValue, double *objectiveFunctionGradient, int *dataIndices, int numDataIndices);

    /**
     * @brief This function is called after each iteration. See IpOpt for arguments.
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
    virtual int intermediateFunction(int alg_mod,
                             int iter_count,
                             double obj_value,
                             double inf_pr, double inf_du,
                             double mu,
                             double d_norm,
                             double regularization_size,
                             double alpha_du, double alpha_pr,
                             int ls_trials);

    /**
     * @brief Called after each cost function evaluation for logging results.
     * @param parameters
     * @param objectiveFunctionValue
     * @param objectiveFunctionGradient
     * @param numFunctionCalls
     * @param timeElapsed
     */
    virtual void logObjectiveFunctionEvaluation(const double *parameters,
                                                double objectiveFunctionValue,
                                                const double *objectiveFunctionGradient,
                                                int numFunctionCalls,
                                                double timeElapsed);

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
                                int exitStatus);

    /**
     * @brief Is called by worker processes to run a simulation for the given condition
     * @param udata UserData for simulation. Model dimensions, sensitivity options and UserData::p is set, others are not.
     * E.g. UserData::k has to be update if applicable.
     * @param dataProvider
     * @param path
     * @param jobId
     * @param resultWriter
     * @param status
     * @return
     */
    static ReturnData *runAndLogSimulation(UserData *udata, MultiConditionDataProvider *dataProvider,
                                           JobIdentifier path, int jobId, MultiConditionProblemResultWriter *resultWriter, int *status);

    virtual MultiConditionDataProvider *getDataProvider();

    virtual double *getInitialParameters(int multiStartIndex) const;

    ~MultiConditionProblem();


    JobIdentifier path;

protected:

    void init();

    void updateUserDataCommon(const double *simulationParameters, const double *objectiveFunctionGradient);

    /**
     * @brief Run AMICI simulations for conditions with the given indices
     * @param optimizationVariables
     * @param logLikelihood
     * @param objectiveFunctionGradient
     * @param dataIndices
     * @param numDataIndices
     * @return Simulation status, != 0 indicates failure
     */
    virtual int runSimulations(const double *optimizationVariables, double *logLikelihood, double *objectiveFunctionGradient, int *dataIndices, int numDataIndices);

    /**
     * @brief Aggregates loglikelihood received from workers.
     * @param data
     * @param logLikelihood
     * @param objectiveFunctionGradient
     * @param dataIndices
     * @param numDataIndices
     * @return *Negative* log likelihood.
     */

    int aggregateLikelihood(JobData *data, double *logLikelihood, double *objectiveFunctionGradient, int *dataIndices, int numDataIndices);

    void printObjectiveFunctionFailureMessage();

    /**
     * @brief Aggregates loglikelihood gradient received from workers.
     * @param conditionIdx
     * @param simulationGradient
     * @param objectiveFunctionGradient
     * @param Gradient of the *negative* log likelihood.
     */

    void addSimulationGradientToObjectiveFunctionGradient(int conditionIdx, const double *simulationGradient, double *objectiveFunctionGradient, int numCommon);

    void addSimulationGradientToObjectiveFunctionGradientConditionSpecificParameters(const double *simulationGradient, double *objectiveFunctionGradient, int numCommon, int numConditionSpecificParams, int firstIndexOfCurrentConditionsSpecificOptimizationParameters);

    int unpackSimulationResult(JobData *d, double *sllhBuffer, double *llh);

    void queueSimulation(JobIdentifier path, JobData *d, int *jobDone,
                         pthread_cond_t *jobDoneChangedCondition, pthread_mutex_t *jobDoneChangedMutex,
                         int lenSendBuffer);

    virtual void setSensitivityOptions(bool sensiRequired);


    UserData *udata;

    MultiConditionDataProvider *dataProvider;

    // keep track of previous results to avoid re-evaluation at the same parameters (using IpOpt new_x)
    double *lastOptimizationParameters;
    double *lastObjectiveFunctionGradient;
    double lastObjectiveFunctionValue;
};



/**
 * @brief The MultiConditionProblemSerial class does the same as its base class, but runs simulations
 * by itself within the same thread. Mostly intended for debugging.
 */

class MultiConditionProblemSerial : public MultiConditionProblem {

public:
    MultiConditionProblemSerial() {}

    MultiConditionProblemSerial(MultiConditionDataProvider *dataProvider) : MultiConditionProblem(dataProvider) {}

    int runSimulations(const double *optimizationVariables, double *logLikelihood, double *objectiveFunctionGradient,  int *dataIndices, int numDataIndices);

};


/**
 * @brief The MultiConditionProblemGeneratorForMultiStart class generates new MultiConditionProblem instances
 * with proper DataProviders for multi-start optimization
 */

class MultiConditionProblemGeneratorForMultiStart : public OptimizationProblemGeneratorForMultiStart {
public:
    OptimizationProblem *getLocalProblemImpl(int multiStartIndex);

    MultiConditionDataProvider *dp;
    OptimizationOptions *options;
    MultiConditionProblemResultWriter *resultWriter;

};

/**
 * @brief Callback function for loadbalancer
 * @param buffer In/out: message buffer
 * @param msgSize In/out: size (bytes) of buffer
 * @param jobId: In: Identifier of the job (unique up to INT_MAX)
 * @param userData In: Pointer to MultiConditionProblem
 */
void handleWorkPackage(char **buffer, int *msgSize, int jobId, void *userData);

#endif
