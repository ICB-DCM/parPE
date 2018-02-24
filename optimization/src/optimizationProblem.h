#ifndef OPTIMIZATION_PROBLEM_H
#define OPTIMIZATION_PROBLEM_H

#include "optimizationResultWriter.h"
#include <optimizationOptions.h>
#include <misc.h>
#include <functions.h>

#include <cstdlib>
#include <vector>
#include <ctime>

#include <hdf5.h>

namespace parpe {

class OptimizationResultWriter;
class OptimizationReporter;


template<typename T>
class SummedGradientProblem {

public:
    SummedGradientProblem() = default;
    SummedGradientProblem(std::unique_ptr<SummedGradientFunction<T>> costFun);

    virtual ~SummedGradientProblem() = default;

    /** Default implementation: random starting points are drawn from [parametersMin, parametersMax] */
    virtual void fillInitialParameters(double *buffer) const;

    /** lower bound of parameter values */
    virtual void fillParametersMin(double *buffer) const = 0;

    /** upper bound of parameter values */
    // TODO:     template <class RandomAccessIterator>
    virtual void fillParametersMax(double *buffer) const = 0;

    OptimizationOptions const& getOptimizationOptions() const;

    void setOptimizationOptions(OptimizationOptions const& options);

    // const?
    std::unique_ptr<SummedGradientFunction<T>> costFun;

    virtual std::unique_ptr<OptimizationReporter> getReporter() const;

private:
    OptimizationOptions optimizationOptions;
};





/**
 * @brief The OptimizationReporter class is called from the optimizer and takes care of
 * things like of logging intermediate results, timing and can tell the optimizer to exit
 */

class OptimizationReporter {
public:
    OptimizationReporter();

    OptimizationReporter(std::unique_ptr<OptimizationResultWriter> rw);

    /**
     * @brief Is called just before the optimizer starts. Must be called before other functions.
     * @param numParameters
     * @param initialParameters
     * @return Quit optimization?
     */
    virtual bool starting(int numParameters, double const *const initialParameters);


    /**
     * @brief Is called after each iteration except for the last one
     * @param numParameters
     * @param parameters
     * @param currentIter
     * @return Quit optimization?
     */
    virtual bool iterationFinished(int numParameters, double const *const parameters, double objectiveFunctionValue,
                                   double const *const objectiveFunctionGradient);

    virtual bool beforeCostFunctionCall(int numParameters, double const *const parameters);

    virtual bool afterCostFunctionCall(int numParameters, double const *const parameters,
                                       double objectiveFunctionValue,
                                       double const *const objectiveFunctionGradient);

    /**
     * @brief Is called after optimization finished
     */
    virtual void finished(double optimalCost,
                          const double *optimalParameters, int exitStatus);


    // TODO how to pass optimizer-specific info? pass OptimizerStatus class ?

    //    virtual int intermediateFunction(int alg_mod, int iter_count,
    //                                     double obj_value, double inf_pr,
    //                                     double inf_du, double mu, double d_norm,
    //                                     double regularization_size,
    //                                     double alpha_du, double alpha_pr,
    //                                     int ls_trials);

private:
    WallTimer wallTimer;
//    clock_t timeOptimizationBegin;
//    clock_t timeIterationBegin;
//    clock_t timeCostEvaluationBegin;

    std::unique_ptr<OptimizationResultWriter> resultWriter;
    int numFunctionCalls = 0;
    int numIterations = 0;
    int numParameters = 0;

    bool started = false;

};



/**
 * @brief The OptimizationProblem class describes an optimization problem.
 *
 * A OptimizationProblem has a GradientFunction objective function to be minimized,
 * parameter bounds and initial values.
 *
 * Additional constraints are currently not supported.
 *
 * TODO: rename GradientProblem? Turn into interface
 * should not have state; so cannot track iterations in there
 */

class OptimizationProblem {

public:
    OptimizationProblem() = default;
    OptimizationProblem(std::unique_ptr<GradientFunction> costFun);

    virtual ~OptimizationProblem() = default;

    /** Default implementation: random starting points are drawn from [parametersMin, parametersMax] */
    virtual void fillInitialParameters(double *buffer) const;

    /** lower bound of parameter values */
    virtual void fillParametersMin(double *buffer) const = 0;

    /** upper bound of parameter values */
    // TODO:     template <class RandomAccessIterator>
    virtual void fillParametersMax(double *buffer) const = 0;

    virtual OptimizationOptions const& getOptimizationOptions() const;

    virtual void setOptimizationOptions(OptimizationOptions const& options);

    // const?
    std::unique_ptr<GradientFunction> costFun;

    virtual std::unique_ptr<OptimizationReporter> getReporter() const;

private:
    OptimizationOptions optimizationOptions;
};


/**
 * @brief Mixin class for handling parameter bounds
 */
class OptimizationProblemImpl : public OptimizationProblem {

public:
    using OptimizationProblem::OptimizationProblem;

    /** lower bound of parameter values */
    void fillParametersMin(double *buffer) const {
        std::copy(parametersMin.begin(), parametersMin.end(), buffer);
    }

    /** upper bound of parameter values */
    void fillParametersMax(double *buffer) const {
        std::copy(parametersMax.begin(), parametersMax.end(), buffer);
    }

    void setParametersMin(std::vector<double> parametersMin) {
        this->parametersMin = parametersMin;
    }

    void setParametersMax(std::vector<double> parametersMax) {
        this->parametersMax = parametersMax;
    }

    void setInitialParameters(std::vector<double> initial) {
        parametersStart = initial;
    }

    void fillInitialParameters(double *buffer) const override {
        if(parametersStart.size()) {
            std::copy(parametersStart.begin(), parametersStart.end(), buffer);
        } else {
            OptimizationProblem::fillInitialParameters(buffer);
        }
    }

private:
    std::vector<double> parametersMin;
    std::vector<double> parametersMax;
    std::vector<double> parametersStart;

};



int getLocalOptimum(OptimizationProblem *problem);

void *getLocalOptimumThreadWrapper(void *optimizationProblemVp);

void runOptimizationsParallel(const OptimizationProblem **problems,
                              int numProblems);

void optimizationProblemGradientCheck(OptimizationProblem *problem,
                                      int numParameterIndicesToCheck, double epsilon);

void optimizationProblemGradientCheck(OptimizationProblem *problem,
                                      const int parameterIndices[],
                                      int numParameterIndices, double epsilon);

} // namespace parpe

#endif
