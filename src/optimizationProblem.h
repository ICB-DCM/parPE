#ifndef OPTIMIZATION_PROBLEM_H
#define OPTIMIZATION_PROBLEM_H

#include <stdbool.h>

typedef enum optimizer_tag {OPTIMIZER_IPOPT, OPTIMIZER_CERES} optimizerEnum;

typedef struct OptimizationProblem_tag OptimizationProblem;

typedef int (*objectiveFunctionFp)(OptimizationProblem *problem,
                                   const double *parameters,
                                   double *result);

typedef int (*objectiveFunctionGradientFp)(OptimizationProblem *problem,
                                            const double *parameters,
                                            double *objFunVal,
                                            double *objFunGrad);

typedef int (*intermediateFunctionFp)(OptimizationProblem *problem,
                                       int alg_mod,
                                       int iter_count,
                                       double obj_value,
                                       double inf_pr, double inf_du,
                                       double mu,
                                       double d_norm,
                                       double regularization_size,
                                       double alpha_du, double alpha_pr,
                                       int ls_trials);

typedef void (*logObjectiveFunctionEvaluationFp)(OptimizationProblem *problem,
                                                 const double *parameters,
                                                 double objectiveFunctionValue,
                                                 int numFunctionCalls,
                                                 double timeElapsed);

typedef void (*logObjectiveFunctionGradientEvaluationFp)(OptimizationProblem *problem,
                                                         const double *parameters,
                                                         double objectiveFunctionValue,
                                                         const double *objectiveFunctionGradient,
                                                         int numFunctionCalls,
                                                         double timeElapsed);
// typedef void (*logOptimizerIterationFp)(void *);

typedef void (*logOptimizerFinishedFp)(OptimizationProblem *problem,
                                       double optimalCost,
                                       const double *optimalParameters,
                                       double masterTime,
                                       int exitStatus);

typedef struct OptimizationProblem_tag {
    objectiveFunctionFp objectiveFunction;
    objectiveFunctionGradientFp objectiveFunctionGradient;

    int numOptimizationParameters;
    double *initialParameters;
    double *parametersMin;
    double *parametersMax;
    void *userData;

    optimizerEnum optimizer;

    char *logFile;
    bool *printToStdout;
    int maxOptimizerIterations;

    intermediateFunctionFp intermediateFunction;
    logObjectiveFunctionEvaluationFp logObjectiveFunctionEvaluation;
    logObjectiveFunctionGradientEvaluationFp logObjectiveFunctionGradientEvaluation;
    //logOptimizerIterationFp logOptimizerIteration;
    logOptimizerFinishedFp logOptimizerFinished;

} OptimizationProblem;

OptimizationProblem *optimizationProblemNew();

int getLocalOptimum(OptimizationProblem *problem);

void *getLocalOptimumThreadWrapper(void *optimizationProblemVp);

#endif
