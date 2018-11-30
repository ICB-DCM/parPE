#ifndef OPTIMIZATION_PROBLEM_H
#define OPTIMIZATION_PROBLEM_H

#include "optimizationResultWriter.h"
#include <optimizationOptions.h>
#include <misc.h>
#include <functions.h>
#include <logging.h>

#include <cstdlib>
#include <vector>
#include <ctime>
#include <cmath>

#include <hdf5.h>
#include <gsl/gsl-lite.hpp>

namespace parpe {

class OptimizationResultWriter;
class OptimizationReporter;

///**
// * TODO: will we use this for minibatch?
// */
//template<typename T>
//class SummedGradientProblem {

//public:
//    SummedGradientProblem() = default;
//    SummedGradientProblem(std::unique_ptr<SummedGradientFunction<T>> costFun);

//    virtual ~SummedGradientProblem() = default;

//    /** Default implementation: random starting points are drawn from [parametersMin, parametersMax] */
//    virtual void fillInitialParameters(double *buffer) const;

//    /** lower bound of parameter values */
//    virtual void fillParametersMin(double *buffer) const = 0;

//    /** upper bound of parameter values */
//    // TODO:     template <class RandomAccessIterator>
//    virtual void fillParametersMax(double *buffer) const = 0;

//    OptimizationOptions const& getOptimizationOptions() const;

//    void setOptimizationOptions(OptimizationOptions const& options);

//    // const?
//    std::unique_ptr<SummedGradientFunction<T>> costFun;

//    virtual std::unique_ptr<OptimizationReporter> getReporter() const;

//private:
//    OptimizationOptions optimizationOptions;
//};


/**
 * @brief The OptimizationReporter class is called from the optimizer and takes care of
 * calling the objective function and things like of keeping track of iterations, logging intermediate results, timing
 * and can tell the optimizer to exit.
 *
 * This extra level of abstraction is added to avoid reimplementing timing, and other things for
 * each supported optimizer.
 */

class OptimizationReporter : public GradientFunction {
public:
    OptimizationReporter(GradientFunction *gradFun, std::unique_ptr<Logger> logger);

    OptimizationReporter(GradientFunction *gradFun,
                         std::unique_ptr<OptimizationResultWriter> rw,
                         std::unique_ptr<Logger> logger);

    virtual ~OptimizationReporter() override = default;

    FunctionEvaluationStatus evaluate(
            gsl::span<double const> parameters,
            double &fval,
            gsl::span<double> gradient,
            Logger *logger = nullptr,
            double *cpuTime = nullptr) const override;

    int numParameters() const override;

    /**
     * @brief Is called just before the optimizer starts. Must be called before other functions.
     * @param numParameters
     * @param initialParameters
     * @return Quit optimization?
     */
    virtual bool starting(gsl::span<const double> initialParameters) const;


    /**
     * @brief Is called after each iteration except for the last one
     * @param numParameters
     * @param parameters
     * @param currentIter
     * @return Quit optimization?
     */
    virtual bool iterationFinished(gsl::span<const double> parameters,
                                   double objectiveFunctionValue,
                                   gsl::span<const double> objectiveFunctionGradient) const;

    virtual bool beforeCostFunctionCall(gsl::span<const double> parameters) const;

    virtual bool afterCostFunctionCall(gsl::span<const double> parameters,
                                       double objectiveFunctionValue,
                                       gsl::span<double const> objectiveFunctionGradient) const;

    /**
     * @brief Is called after optimization finished
     */
    virtual void finished(double optimalCost,
                          gsl::span<const double> parameters, int exitStatus) const;


    // TODO how to pass optimizer-specific info? pass OptimizerStatus class ?

    //    virtual int intermediateFunction(int alg_mod, int iter_count,
    //                                     double obj_value, double inf_pr,
    //                                     double inf_du, double mu, double d_norm,
    //                                     double regularization_size,
    //                                     double alpha_du, double alpha_pr,
    //                                     int ls_trials);

    virtual double getFinalCost() const;

    virtual std::vector<double> const& getFinalParameters() const;

    void setGradientFunction(GradientFunction *gradFun) const;

    void setLoggingEachIteration(bool logGradient) const;
    
    void setLoggingEachFunctionEvaluation(bool logGradient, bool logParameters) const;
    
    std::unique_ptr<OptimizationResultWriter> resultWriter;

    mutable double cpuTimeTotalSec = 0.0;
    mutable double cpuTimeIterationSec = 0.0;
    std::unique_ptr<Logger> logger;

protected:
    void printObjectiveFunctionFailureMessage() const;


    // data members are mutable, because we inherit from GradientFunction,
    // and evaluate() is const there. This could probably be solved better....

    mutable WallTimer wallTimer;

//    clock_t timeOptimizationBegin;
//    clock_t timeIterationBegin;
//    clock_t timeCostEvaluationBegin;

    mutable int numFunctionCalls = 0;
    mutable int numIterations = 0;
    mutable int numParameters_ = 0;

    mutable bool started = false;
    
    mutable bool logGradientEachIteration = true;
    mutable bool logGradientEachFunctionEvaluation = true;
    mutable bool logParametersEachFunctionEvaluation = true;

    // non-owning
    mutable GradientFunction *gradFun = nullptr;

    // for caching
    mutable bool haveCachedCost = false;
    mutable bool haveCachedGradient = false;
    mutable std::vector<double> cachedGradient;
    mutable double cachedCost = std::numeric_limits<double>::infinity();
    mutable FunctionEvaluationStatus cachedStatus = functionEvaluationSuccess;

    mutable double finalCost = std::numeric_limits<double>::quiet_NaN();

    // keeps the most recent parameters, assuming they are the final ones
    mutable std::vector<double> cachedParameters;

    std::string defaultLoggerPrefix;
private:

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
 */

class OptimizationProblem {

public:
    OptimizationProblem() = default;
    OptimizationProblem(std::unique_ptr<GradientFunction> costFun,
                        std::unique_ptr<Logger> logger);
    OptimizationProblem(OptimizationProblem const& other) = delete;

    virtual ~OptimizationProblem() = default;

    /** Default implementation: random starting points are drawn from [parametersMin, parametersMax] */
    virtual void fillInitialParameters(gsl::span<double> buffer) const;

    /** lower bound of parameter values */
    virtual void fillParametersMin(gsl::span<double> buffer) const = 0;

    /** upper bound of parameter values */
    // TODO:     template <class RandomAccessIterator>
    virtual void fillParametersMax(gsl::span<double> buffer) const = 0;

    virtual OptimizationOptions const& getOptimizationOptions() const;

    virtual void setOptimizationOptions(OptimizationOptions const& options);

    virtual std::unique_ptr<OptimizationReporter> getReporter() const;

    // const?
    std::unique_ptr<GradientFunction> costFun;

    std::unique_ptr<Logger> logger;

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
    void fillParametersMin(gsl::span<double> buffer) const override {
        std::copy(parametersMin.begin(), parametersMin.end(), buffer.begin());
    }

    /** upper bound of parameter values */
    void fillParametersMax(gsl::span<double> buffer) const override {
        std::copy(parametersMax.begin(), parametersMax.end(), buffer.begin());
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

    void fillInitialParameters(gsl::span<double> buffer) const override {
        if(parametersStart.size()) {
            std::copy(parametersStart.begin(), parametersStart.end(), buffer.begin());
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

//void runOptimizationsParallel(OptimizationProblem **problems,
//                              int numProblems);

void optimizationProblemGradientCheck(OptimizationProblem *problem,
                                      int numParameterIndicesToCheck, double epsilon);

void optimizationProblemGradientCheck(OptimizationProblem *problem,
                                      gsl::span<const int> parameterIndices,
                                      double epsilon);

} // namespace parpe

#endif
