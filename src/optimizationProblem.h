#ifndef OPTIMIZATION_PROBLEM_H
#define OPTIMIZATION_PROBLEM_H

#include <stdbool.h>

typedef enum optimizer_tag {OPTIMIZER_IPOPT, OPTIMIZER_CERES} optimizerEnum;

typedef struct OptimizationProblem_tag OptimizationProblem;

/**
 * Callback function for objective function evaluation at parameters.
 * Non-zero return status indicates failure.
 */
typedef int (*objectiveFunctionFp)(OptimizationProblem *problem,
                                   const double *parameters,
                                   double *result);

/**
 * Callback function for objective function gradient evaluation at parameters.
 * Non-zero return status indicates failure.
 */
typedef int (*objectiveFunctionGradientFp)(OptimizationProblem *problem,
                                            const double *parameters,
                                            double *objFunVal,
                                            double *objFunGrad);

/**
 * Callback function which is called after each optimizer iteration.
 * Non-zero return status indicates failure.
 */
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

// TODO: remove callback functions; let user handle within objective function evaluation
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

typedef void (*logOptimizerFinishedFp)(OptimizationProblem *problem,
                                       double optimalCost,
                                       const double *optimalParameters,
                                       double masterTime,
                                       int exitStatus);

/** Type to describe an optimization (minimization) problem */

typedef struct OptimizationProblem_tag {
    /** pointer to cost function */
    objectiveFunctionFp objectiveFunction;
    /** pointer to gradient function of the cost function */
    objectiveFunctionGradientFp objectiveFunctionGradient;

    /** number of optimization parameters */
    int numOptimizationParameters;
    /** starting point for optimization */
    double *initialParameters;
    /** lowest allowed parameter values */
    double *parametersMin;
    /** highest allowed parameter values */
    double *parametersMax;
    /** any user-provided data that is passed to all callback functions */
    void *userData;

    /** Optimizer to use */
    optimizerEnum optimizer;

    /** Optimizer log file */
    char *logFile;

    /** Print progress to stdout */
    bool printToStdout;
    /** Maximum number of optimizer iterations*/
    int maxOptimizerIterations;

    /** function to be called after each optimizer iteration */
    intermediateFunctionFp intermediateFunction;

    // TODO remove
    logObjectiveFunctionEvaluationFp logObjectiveFunctionEvaluation;
    logObjectiveFunctionGradientEvaluationFp logObjectiveFunctionGradientEvaluation;
    logOptimizerFinishedFp logOptimizerFinished;

} OptimizationProblem;

/** Create a new optimization problem with default values */
OptimizationProblem *optimizationProblemNew();

int getLocalOptimum(OptimizationProblem *problem);

void *getLocalOptimumThreadWrapper(void *optimizationProblemVp);

void runOptimizationsParallel(const OptimizationProblem **problems, int numProblems);
#endif
