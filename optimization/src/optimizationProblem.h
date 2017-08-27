#ifndef OPTIMIZATION_PROBLEM_H
#define OPTIMIZATION_PROBLEM_H

#include "optimizationResultWriter.h"
#include <cstdlib>
#include <hdf5.h>

class OptimizationOptions;
class OptimizationResultWriter;

/**
 * @brief The OptimizationProblem class describes an optimization problem
 */
class OptimizationProblem {

  public:
    /**
     * Callback function for objective function gradient evaluation at
     * parameters.
     * Non-zero return status indicates failure. If objFunGrad is not null,
     * gradient information is expected.
     */
    virtual int evaluateObjectiveFunction(const double *parameters,
                                          double *objFunVal,
                                          double *objFunGrad) = 0;

    /**
     * Callback function which is called after each optimizer iteration.
     * Non-zero return status indicates failure.
     */
    virtual int intermediateFunction(int alg_mod, int iter_count,
                                     double obj_value, double inf_pr,
                                     double inf_du, double mu, double d_norm,
                                     double regularization_size,
                                     double alpha_du, double alpha_pr,
                                     int ls_trials);

    virtual void
    logObjectiveFunctionEvaluation(const double *parameters,
                                   double objectiveFunctionValue,
                                   const double *objectiveFunctionGradient,
                                   int numFunctionCalls, double timeElapsed);

    virtual void logOptimizerFinished(double optimalCost,
                                      const double *optimalParameters,
                                      double masterTime, int exitStatus);

    virtual ~OptimizationProblem();

    virtual double *getInitialParameters(int multiStartIndex) const;

    virtual double *getInitialParameters() const;

    /** random starting points are drawn from [parametersMin, parametersMax] */

    void fillInitialParameters(double *buffer) const;

    static void getRandomStartingpoint(const double *min, const double *max,
                                       int numParameters, double *buffer);

    void setInitialParameters(double *initialParameters);

    int getNumOptimizationParameters() const;

    const double *getParametersMin() const;

    const double *getParametersMax() const;

    OptimizationOptions *optimizationOptions = nullptr;

  protected:
    /** number of optimization parameters */
    int numOptimizationParameters = 0;

    /** lowest allowed parameter values */
    double *parametersMin = nullptr;

    /** highest allowed parameter values */
    double *parametersMax = nullptr;

    double *initialParameters = nullptr;
};

int getLocalOptimum(OptimizationProblem *problem);

void *getLocalOptimumThreadWrapper(void *optimizationProblemVp);

void runOptimizationsParallel(const OptimizationProblem **problems,
                              int numProblems);

void optimizationProblemGradientCheck(OptimizationProblem *problem,
                                      const int parameterIndices[],
                                      int numParameterIndices, double epsilon);

#endif
