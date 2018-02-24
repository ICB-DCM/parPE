#include <CppUTest/TestHarness.h>
#include <CppUTestExt/MockSupport.h>

#include "optimizationProblem.h"
#include "quadraticTestProblem.h"
#include <misc.h>
#include "localOptimizationIpopt.h"

#include <cmath>

// clang-format off
TEST_GROUP(optimizationProblem){
    void setup() {
        mock().disable();
    }

    void teardown() {
    }
};
// clang-format on


/**
 * @brief The SummedGradientFunctionLinearModelTest class is a linear model with
 * mean squared error cost function.
 *
 * y = a * x + b
 *
 * cost = MSE = 1/N \sum_i^N (\bar{y} - y)^2
 */
class SummedGradientFunctionLinearModelTest : public parpe::SummedGradientFunction<double> {
public:
    virtual parpe::FunctionEvaluationStatus evaluate(
            const double* const parameters,
            double dataset,
            double &fval,
            double* gradient) const override
    {
        fval = parameters[0] * dataset + parameters[1];

        if(gradient) {
            gradient[0] = dataset;
            gradient[1] = 0;
        }

        return parpe::functionEvaluationSuccess;
    }

    /**
     * @brief Evaluate cost function on a set of training points. Calls evaluate for every
     * data points. This is far from efficient and intended to be used only for testing.
     * @param parameters
     * @param datasets
     * @param fval
     * @param gradient
     * @return
     */
    virtual parpe::FunctionEvaluationStatus evaluate(
            const double* const parameters,
            std::vector<double> datasets,
            double &fval,
            double* gradient) const override
    {
        fval = 0;
        if(gradient)
            std::fill(gradient, gradient + numParameters(), 0);

        double tmpFVal = 0;
        std::vector<double> tmpGradient(numParameters());

        for(auto& d : datasets) {
            auto status = evaluate(parameters, d, tmpFVal, tmpGradient.data());
            if(status != parpe::functionEvaluationSuccess)
                return status;

            fval += tmpFVal;
            if(gradient) {
                for(int i = 0; i < numParameters(); ++i)
                    gradient[i] += tmpGradient[i];
            }
        }
        return parpe::functionEvaluationSuccess;
    }

    virtual int numParameters() const override
    {
        return numParameters_;
    }

    int numParameters_ = 2;

};


TEST(optimizationProblem, quadraticTestFunction) {
    // Test QuadraticGradientFunction for f(-1) = 42
    parpe::QuadraticGradientFunction f {};
    double parameter = -1;
    double fValExp = 42.0;
    double gradientExp = 0;

    double fValAct = NAN;
    double gradientAct = NAN;
    f.evaluate(&parameter, fValAct, &gradientAct);

    CHECK_EQUAL(fValExp, fValAct);
    CHECK_EQUAL(gradientExp, gradientAct);
}



TEST(optimizationProblem, gradientChecker) {
    parpe::QuadraticTestProblem problem {};
    constexpr int numParameterIndices {1};
    int parameterIndices[numParameterIndices] {0};

    parpe::optimizationProblemGradientCheck(&problem, parameterIndices, numParameterIndices, 1e-6);
}


TEST(optimizationProblem, linearModel) {
    // Test if the linear model produces correct results
    SummedGradientFunctionLinearModelTest model;

    std::vector<double> parameters {1.0, 2.0};
    double fval = NAN;
    std::vector<double> gradient(model.numParameters(), NAN);

    model.evaluate(parameters.data(), 1.0, fval, gradient.data());
    CHECK_EQUAL(3.0, fval);
    CHECK_EQUAL(1.0, gradient[0]);
    CHECK_EQUAL(0.0, gradient[1]);

    std::vector<double> dataset {2.0, 3.0};
    model.evaluate(parameters.data(), dataset, fval, gradient.data());
    CHECK_EQUAL(9.0, fval);
    CHECK_EQUAL(5.0, gradient[0]);
    CHECK_EQUAL(0.0, gradient[1]);
}


TEST(optimizationProblem, linearModelToGradientFun) {
    // Test that the SummedGradientFunction <-> GradientFunction Adapter works with the linear model
    std::vector<double> dataset {2.0, 3.0};
    auto model = std::unique_ptr<parpe::SummedGradientFunction<double>>(new SummedGradientFunctionLinearModelTest());
    parpe::SummedGradientFunctionGradientFunctionAdapter<double> gradFun(std::move(model), dataset);

    std::vector<double> parameters {1.0, 2.0};
    double fval = NAN;
    std::vector<double> gradient(gradFun.numParameters(), NAN);

    gradFun.evaluate(parameters.data(), fval, gradient.data());
    CHECK_EQUAL(9.0, fval);
    CHECK_EQUAL(5.0, gradient[0]);
    CHECK_EQUAL(0.0, gradient[1]);
}



TEST(optimizationProblem, linearModelToGradientFunOptimization) {
    // create optimization problem for the linear model
    // does not do anything meaningful yet
    std::vector<double> dataset {2.0, 3.0};
    auto model = std::unique_ptr<parpe::SummedGradientFunction<double>>(new SummedGradientFunctionLinearModelTest());
    auto gradFun = std::unique_ptr<parpe::GradientFunction>(
                new parpe::SummedGradientFunctionGradientFunctionAdapter<double>(
                    std::move(model), dataset));
    parpe::OptimizationProblemImpl problem {std::move(gradFun)};

    std::vector<double> parametersMin {-100, -100};
    std::vector<double> parametersMax {100, 100};

    problem.setParametersMin(parametersMin);
    problem.setParametersMax(parametersMax);

    parpe::OptimizerIpOpt opt;
    opt.optimize(&problem);
}
