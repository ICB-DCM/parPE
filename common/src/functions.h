#ifndef PARPE_OPTIMIZATION_FUNCTIONS_H
#define PARPE_OPTIMIZATION_FUNCTIONS_H

#include <logging.h>

#include <vector>
#include <memory>

#include <gsl/gsl-lite.hpp>

namespace parpe {



enum FunctionEvaluationStatus {
    functionEvaluationSuccess,
    functionEvaluationFailure,
};


/**
 * @brief The GradientFunction class is an interface for an
 * arbitrary function f(x) and its gradient.
 */

class GradientFunction {
public:

    /**
     * @brief Evaluate the function f(x)
     * @param parameters Point x at which to evaluate f(x). Must be of length numParameters().
     * @param fval (output) Will be set to the function value f(x)
     * @param gradient (output) If not nullptr, will contain the gradient of f(x) at x. Must be of length numParameters().
     * @return functionEvaluationSuccess on success, functionEvaluationFailure otherwise
     */
    virtual FunctionEvaluationStatus evaluate(
            gsl::span<double const> parameters,
            double &fval,
            gsl::span<double> gradient) const;


    virtual FunctionEvaluationStatus evaluate(
            gsl::span<double const> parameters,
            double &fval,
            gsl::span<double> gradient,
            Logger *logger,
            double *cpuTime) const = 0;

    virtual int numParameters() const = 0;

    virtual ~GradientFunction() = default;
};



/**
 * @brief The SummedGradientFunction class is an interface for cost functions and gradients that are a sum of functions evaluated on a number of data records.
 * To be used e.g. for mini-batch optimization. Template parameter can be used for data indices
 * or directly for data points.
 */
template<typename T>
class SummedGradientFunction {
public:
    /**
     * @brief Evaluate on single data point
     * @param parameters
     * @param dataset
     * @param fval
     * @param gradient
     * @return
     */
    virtual FunctionEvaluationStatus evaluate(
            gsl::span<const double> parameters,
            T dataset,
            double &fval,
            gsl::span<double> gradient,
            Logger *logger,
            double *cpuTime) const = 0;

    /**
     * @brief Evaluate on vector of data points
     * @param parameters
     * @param datasets
     * @param fval
     * @param gradient
     * @return
     */
    virtual FunctionEvaluationStatus evaluate(
            gsl::span<const double> parameters,
            std::vector<T> datasets,
            double &fval,
            gsl::span<double> gradient,
            Logger *logger,
            double *cpuTime
            ) const = 0;

    virtual int numParameters() const = 0;

    virtual ~SummedGradientFunction() = default;
};



/**
 * Adapter / wrapper for SummedGradientFunction to GradientFunction..
 *
 * Simply evaluates SummedGradientFunction on all datasets.
 */
template<typename T>
class SummedGradientFunctionGradientFunctionAdapter
        : public GradientFunction, public SummedGradientFunction<T>
{
public:
    /**
     * @brief SummedGradientFunctionGradientFunctionAdapter
     * @param gradFun Function to be wrapped
     * @param datasets Datasets on which to evaluate
     */
    SummedGradientFunctionGradientFunctionAdapter(
            std::unique_ptr< SummedGradientFunction<T> > gradFun,
            std::vector<T> datasets)
        : gradFun(std::move(gradFun)),
          datasets(datasets)
    {
    }

    FunctionEvaluationStatus evaluate(
            gsl::span<const double> parameters,
            double &fval,
            gsl::span<double> gradient,
            Logger *logger = nullptr,
            double *cpuTime = nullptr) const override
    {
        return gradFun->evaluate(parameters, datasets, fval, gradient, logger, cpuTime);
    }

    FunctionEvaluationStatus evaluate(
            gsl::span<const double> parameters,
            T dataset,
            double &fval,
            gsl::span<double> gradient,
            Logger *logger,
            double *cpuTime) const override
    {
        return gradFun->evaluate(parameters, dataset, fval, gradient, logger, cpuTime);
    }

    FunctionEvaluationStatus evaluate(
            gsl::span<const double> parameters,
            std::vector<T> datasets,
            double &fval,
            gsl::span<double> gradient,
            Logger *logger,
            double *cpuTime
            ) const override
    {
        return gradFun->evaluate(parameters, datasets, fval, gradient, logger, cpuTime);
    }

    int numParameters() const override { return gradFun->numParameters(); }

    SummedGradientFunction<T>* getWrappedFunction() const { return gradFun.get(); }

private:
    std::unique_ptr<SummedGradientFunction<T>> gradFun;
    std::vector<T> datasets;
};


}
#endif
