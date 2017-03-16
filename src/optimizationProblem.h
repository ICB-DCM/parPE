#ifndef OPTIMIZATION_PROBLEM_H
#define OPTIMIZATION_PROBLEM_H


typedef int (*objectiveFunctionFp)(void *problem,
                                   const double *parameters,
                                   double *result);

typedef int (*objectiveFunctionGradientFp)(void *problem,
                                            const double *parameters,
                                            double *objFunVal,
                                            double *objFunGrad);

typedef int (*intermediateFunctionFp)(void *problem,
                                       int alg_mod,
                                       int iter_count,
                                       double obj_value,
                                       double inf_pr, double inf_du,
                                       double mu,
                                       double d_norm,
                                       double regularization_size,
                                       double alpha_du, double alpha_pr,
                                       int ls_trials);

typedef void (*logObjectiveFunctionEvaluationFp)(void *problem,
                                                 const double *parameters,
                                                 double objectiveFunctionValue,
                                                 int numFunctionCalls,
                                                 double timeElapsed);

typedef void (*logObjectiveFunctionGradientEvaluationFp)(void *problem,
                                                         const double *parameters,
                                                         double objectiveFunctionValue,
                                                         const double *objectiveFunctionGradient,
                                                         int numFunctionCalls,
                                                         double timeElapsed);
// typedef void (*logOptimizerIterationFp)(void *);

typedef void (*logOptimizerFinishedFp)(void *problem,
                                       double finalObjectiveFunctionValue,
                                       double masterTime,
                                       int exitStatus);

typedef struct OptimizationProblem_tag {
    objectiveFunctionFp objectiveFunction;
    objectiveFunctionGradientFp objectiveFunctionGradient;
    intermediateFunctionFp intermediateFunction;
    int numOptimizationParameters;
    double *initialParameters;
    double *parametersMin;
    double *parametersMax;
    void *userData;
    char *logFile;
    int maxOptimizerIterations;
    logObjectiveFunctionEvaluationFp logObjectiveFunctionEvaluation;
    logObjectiveFunctionGradientEvaluationFp logObjectiveFunctionGradientEvaluation;
    //logOptimizerIterationFp logOptimizerIteration;
    logOptimizerFinishedFp logOptimizerFinished;
} OptimizationProblem;

OptimizationProblem *optimizationProblemNew();

#endif
