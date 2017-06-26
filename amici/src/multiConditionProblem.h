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

    virtual int evaluateObjectiveFunction(const double *optimiziationVariables, double *objFunVal, double *objectiveFunctionGradient);

    virtual int intermediateFunction(int alg_mod,
                             int iter_count,
                             double obj_value,
                             double inf_pr, double inf_du,
                             double mu,
                             double d_norm,
                             double regularization_size,
                             double alpha_du, double alpha_pr,
                             int ls_trials);

    virtual void logObjectiveFunctionEvaluation(const double *parameters,
                                                double objectiveFunctionValue,
                                                const double *objectiveFunctionGradient,
                                                int numFunctionCalls,
                                                double timeElapsed);

    virtual void logOptimizerFinished(double optimalCost,
                                const double *optimalParameters,
                                double masterTime,
                                int exitStatus);


    static ReturnData *runAndLogSimulation(UserData *udata, MultiConditionDataProvider *dataProvider,
                                           JobIdentifier path, int jobId, OptimizationResultWriter *resultWriter, int *status);

    virtual MultiConditionDataProvider *getDataProvider();

    ~MultiConditionProblem();

    JobIdentifier path;

protected:

    void init();

    void updateUserData(const double *simulationParameters, const double *objectiveFunctionGradient);

    virtual int runSimulations(const double *optimizationVariables, double *logLikelihood, double *objectiveFunctionGradient);

    int aggregateLikelihood(JobData *data, double *logLikelihood, double *objectiveFunctionGradient);

    void printObjectiveFunctionFailureMessage();

    void addSimulationGradientToObjectiveFunctionGradient(int conditionIdx, const double *simulationGradient, double *objectiveFunctionGradient);

    int unpackSimulationResult(JobData *d, double *sllhBuffer, double *llh);

    void queueSimulation(JobIdentifier path, JobData *d, int *jobDone,
                         pthread_cond_t *jobDoneChangedCondition, pthread_mutex_t *jobDoneChangedMutex,
                         int lenSendBuffer);

    void updateUserDataConditionSpecificParameters(int celllineIdx, const double *optimizationParams);

    virtual void setSensitivityOptions(bool sensiRequired);

    UserData *udata;
    MultiConditionDataProvider *dataProvider;

    double *lastOptimizationParameters;
    double *lastObjectiveFunctionGradient;
    double lastObjectiveFunctionValue;
};

class MultiConditionProblemSerial : public MultiConditionProblem {

public:
    MultiConditionProblemSerial() {}

    MultiConditionProblemSerial(MultiConditionDataProvider *dataProvider) : MultiConditionProblem(dataProvider) {}

    int runSimulations(const double *optimizationVariables, double *logLikelihood, double *objectiveFunctionGradient);

};

class MultiConditionProblemGeneratorForMultiStart : public OptimizationProblemGeneratorForMultiStart {
public:
    OptimizationProblem *getLocalProblemImpl(int multiStartIndex);

    MultiConditionDataProvider *dp;
    OptimizationOptions *options;
    MultiConditionProblemResultWriter *resultWriter;

};

void handleWorkPackage(char **buffer, int *msgSize, int jobId, void *userData);

#endif
