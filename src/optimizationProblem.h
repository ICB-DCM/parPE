#ifndef OPTIMIZATION_PROBLEM_H
#define OPTIMIZATION_PROBLEM_H

#include <cstdlib>

typedef enum optimizer_tag {OPTIMIZER_IPOPT, OPTIMIZER_CERES} optimizerEnum;

/** Type to describe an optimization (minimization) problem */

class OptimizationProblem {

public:
    /**
     * Callback function for objective function gradient evaluation at parameters.
     * Non-zero return status indicates failure. If objFunGrad is not null, gradient information is expected.
     */
    virtual int evaluateObjectiveFunction(const double *parameters, double *objFunVal, double *objFunGrad);

    /**
     * Callback function which is called after each optimizer iteration.
     * Non-zero return status indicates failure.
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

    virtual void logObjectiveFunctionEvaluation(const double *parameters,
                                                double objectiveFunctionValue,
                                                const double *objectiveFunctionGradient,
                                                int numFunctionCalls,
                                                double timeElapsed);

    virtual void logOptimizerFinished(double optimalCost,
                                const double *optimalParameters,
                                double masterTime,
                                int exitStatus);

    virtual ~OptimizationProblem(){}

    /** number of optimization parameters */
    int numOptimizationParameters = 0;

    /** starting point for optimization. If 0, random starting points are drawn from [parametersMin, parametersMax] */
    double *initialParameters = NULL;

    /** lowest allowed parameter values */
    double *parametersMin = NULL;

    /** highest allowed parameter values */
    double *parametersMax = NULL;

    /** Optimizer to use */
    optimizerEnum optimizer = OPTIMIZER_IPOPT;

    /** Optimizer log file */
    char *logFile = NULL;

    /** Print progress to stdout */
    bool printToStdout = true;

    /** Maximum number of optimizer iterations*/
    int maxOptimizerIterations = 1000;
};

int getLocalOptimum(OptimizationProblem *problem);

void *getLocalOptimumThreadWrapper(void *optimizationProblemVp);

void runOptimizationsParallel(const OptimizationProblem **problems, int numProblems);

void getRandomStartingpoint(const double *min, const double *max, int numParameters, double *buffer);

void optimizationProblemGradientCheck(OptimizationProblem *problem, const int parameterIndices[], int numParameterIndices, double epsilon);

#endif
