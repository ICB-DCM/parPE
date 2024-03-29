#ifndef OPTIMIZATION_PROBLEM_H
#define OPTIMIZATION_PROBLEM_H

#include <parpeoptimization/optimizationOptions.h>
#include <parpeoptimization/optimizationResultWriter.h>
#include <parpecommon/functions.h>
#include <parpecommon/logging.h>

#include <vector>

#include <hdf5.h>
#include <gsl/gsl-lite.hpp>

namespace parpe {

class OptimizationReporter;

/**
 * @brief The OptimizationReporter class is called from the optimizer and takes
 * care of calling the actual objective function, thereby keeping track of
 * iterations, computation time, logging intermediate results, timing and can
 * tell the optimizer to exit.
 *
 * This extra level of abstraction is added to avoid reimplementing timing, and
 * other things for each supported optimizer. The indirection of cost function
 * evaluation is added to allow caching previous cost function values.
 */

class OptimizationReporter: public GradientFunction {
public:
    OptimizationReporter(GradientFunction *gradFun,
                         std::unique_ptr<Logger> logger);

    OptimizationReporter(GradientFunction *gradFun,
                         std::unique_ptr<OptimizationResultWriter> rw,
                         std::unique_ptr<Logger> logger);

    FunctionEvaluationStatus evaluate(gsl::span<double const> parameters,
                                      double &fval,
                                      gsl::span<double> gradient,
                                      Logger *logger = nullptr,
                                      double *cpuTime = nullptr) const override;

    int numParameters() const override;

    /**
     * @brief Is called just before the optimizer starts. Must be called before
     * other functions.
     * @param initialParameters
     * @return Quit optimization?
     */
    virtual bool starting(gsl::span<const double> initialParameters) const;

    /**
     * @brief Is called after each iteration except for the last one
     * @param parameters
     * @param objectiveFunctionValue
     * @param objectiveFunctionGradient
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
                          gsl::span<const double> parameters,
                          int exitStatus) const;

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

    std::vector<std::string> getParameterIds() const override;

    std::unique_ptr<OptimizationResultWriter> result_writer_;

    mutable double cpu_time_total_sec_ = 0.0;
    mutable double cpu_time_iteration_sec_ = 0.0;
    std::unique_ptr<Logger> logger_;

protected:
    void printObjectiveFunctionFailureMessage() const;

    // data members are mutable, because we inherit from GradientFunction,
    // and evaluate() is const there. This could probably be solved better....

    mutable WallTimer wall_timer_;

    mutable int num_function_calls_ = 0;
    mutable int num_iterations_ = 0;
    mutable int num_parameters_ = 0;

    mutable bool started_ = false;

    // non-owning
    mutable GradientFunction *grad_fun_ = nullptr;

    // for caching
    mutable bool have_cached_cost_ = false;
    mutable bool have_cached_gradient_ = false;
    mutable std::vector<double> cached_gradient_;
    mutable double cached_cost_ = std::numeric_limits<double>::infinity();
    mutable FunctionEvaluationStatus cached_status_ = functionEvaluationSuccess;

    mutable double final_cost_ = std::numeric_limits<double>::quiet_NaN();

    // keeps the most recent parameters, assuming they are the final ones
    mutable std::vector<double> cached_parameters_;

    std::string default_logger_prefix_;
};


/**
 * @brief The OptimizationProblem class describes an optimization problem.
 *
 * A OptimizationProblem has a GradientFunction objective function to be
 * minimized, parameter bounds and initial values.
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

    /** Default implementation: random starting points are drawn from
     * [parametersMin, parametersMax] */
    virtual void fillInitialParameters(gsl::span<double> buffer) const;

    /** lower bound of parameter values */
    virtual void fillParametersMin(gsl::span<double> buffer) const = 0;

    /** upper bound of parameter values */
    virtual void fillParametersMax(gsl::span<double> buffer) const = 0;

    virtual OptimizationOptions const& getOptimizationOptions() const;

    virtual void setOptimizationOptions(OptimizationOptions const& options);

    virtual std::unique_ptr<OptimizationReporter> getReporter() const;

    // const?
    std::unique_ptr<GradientFunction> cost_fun_;

    std::unique_ptr<Logger> logger_;

private:
    OptimizationOptions optimization_options_;
};


/**
 * @brief Mixin class for handling parameter bounds
 */
class OptimizationProblemImpl: public OptimizationProblem {

public:
    using OptimizationProblem::OptimizationProblem;

    /** lower bound of parameter values */
    void fillParametersMin(gsl::span<double> buffer) const override;

    /** upper bound of parameter values */
    void fillParametersMax(gsl::span<double> buffer) const override;

    void setParametersMin(std::vector<double> parametersMin);

    void setParametersMax(std::vector<double> parametersMax);

    void setInitialParameters(std::vector<double> initial);

    void fillInitialParameters(gsl::span<double> buffer) const override;

private:
    std::vector<double> parametersMin;
    std::vector<double> parametersMax;
    std::vector<double> parametersStart;

};


/**
 * @brief getLocalOptimum
 * @param problem
 * @return int indicating status. 0: success, != 0: failure
 */

int getLocalOptimum(OptimizationProblem *problem);


void optimizationProblemGradientCheckMultiEps(OptimizationProblem *problem,
                                              int numParameterIndicesToCheck);

void optimizationProblemGradientCheckMultiEps(OptimizationProblem *problem,
                                              gsl::span<const int> parameterIndices,
                                              gsl::span<double> multi_eps);

void optimizationProblemGradientCheck(OptimizationProblem *problem,
                                      int numParameterIndicesToCheck,
                                      double epsilon);

void optimizationProblemGradientCheck(OptimizationProblem *problem,
                                      gsl::span<const int> parameterIndices,
                                      double epsilon);

} // namespace parpe

#endif
