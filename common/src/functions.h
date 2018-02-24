#ifndef PARPE_OPTIMIZATION_FUNCTIONS_H
#define PARPE_OPTIMIZATION_FUNCTIONS_H

#include <vector>
#include <memory>

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
            const double* const parameters,
            double &fval,
            double* gradient) const = 0;

    virtual int numParameters() const = 0;

    virtual ~GradientFunction() = default;
};



/**
 * @brief The SummedGradientFunction class is an interface for cost functions and gradients that are a sum of functions evaluated on a number of data records.
 * To be used e.g. for mini-batch optimization
 */

template<typename T>
class SummedGradientFunction {
public:
    virtual FunctionEvaluationStatus evaluate(
            const double* const parameters,
            T dataset,
            double &fval,
            double* gradient) const = 0;

    // TODO provide default implementation via function above
    virtual FunctionEvaluationStatus evaluate(
            const double* const parameters,
            std::vector<T> datasets,
            double &fval,
            double* gradient) const = 0;

    virtual int numParameters() const = 0;

    virtual ~SummedGradientFunction() = default;
};




template<typename T>
class SummedGradientFunctionGradientFunctionAdapter : public GradientFunction {
public:
    SummedGradientFunctionGradientFunctionAdapter(std::unique_ptr< SummedGradientFunction<T> > gradFun,
                                                  std::vector<T> datasets)
        : gradFun(std::move(gradFun)), datasets(datasets)
    {
    }

    FunctionEvaluationStatus evaluate(
            const double* const parameters,
            double &fval,
            double* gradient) const override {
        return gradFun->evaluate(parameters, datasets, fval, gradient);

    }

    int numParameters() const override { return gradFun->numParameters(); }

    std::unique_ptr<SummedGradientFunction<T>> gradFun;
private:
    std::vector<T> datasets;
};


}
#endif
