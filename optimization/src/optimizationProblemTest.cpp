#include <CppUTest/TestHarness.h>
#include <CppUTestExt/MockSupport.h>

#include "parpeConfig.h"
#include "optimizationProblem.h"
#include "quadraticTestProblem.h"
#include <misc.h>

#ifdef PARPE_ENABLE_IPOPT
#include "localOptimizationIpopt.h"
#endif

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
    parpe::FunctionEvaluationStatus evaluate(
            gsl::span<double const> parameters,
            double dataset,
            double &fval,
            gsl::span<double> gradient, parpe::Logger* logger, double *cpuTime) const override
    {
        fval = parameters[0] * dataset + parameters[1];

        if(!gradient.empty()) {
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
    parpe::FunctionEvaluationStatus evaluate(
            gsl::span<double const> parameters,
            std::vector<double> datasets,
            double &fval,
            gsl::span<double> gradient,
            parpe::Logger* logger, double* cpuTime) const override
    {
        fval = 0;
        if(!gradient.empty())
            std::fill(gradient.begin(), gradient.end(), 0);

        double tmpFVal = 0;
        std::vector<double> tmpGradient(parameters.size());

        for(auto& d : datasets) {
            auto status = evaluate(parameters, d, tmpFVal, tmpGradient, nullptr, nullptr);
            if(status != parpe::functionEvaluationSuccess)
                return status;

            fval += tmpFVal;
            if(!gradient.empty()) {
                for(int i = 0; i < numParameters(); ++i)
                    gradient[i] += tmpGradient[i];
            }
        }
        return parpe::functionEvaluationSuccess;
    }

    int numParameters() const override
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
    f.evaluate(gsl::span<double const>(&parameter, 1), fValAct, gsl::span<double>(&gradientAct, 1));

    CHECK_EQUAL(fValExp, fValAct);
    CHECK_EQUAL(gradientExp, gradientAct);
}



TEST(optimizationProblem, gradientChecker) {
    parpe::QuadraticTestProblem problem {};
    constexpr int numParameterIndices {1};
    int parameterIndices[numParameterIndices] {0};

    parpe::optimizationProblemGradientCheck(&problem, parameterIndices, 1e-6);
}


TEST(optimizationProblem, linearModel) {
    // Test if the linear model produces correct results
    SummedGradientFunctionLinearModelTest model;

    std::vector<double> parameters {1.0, 2.0};
    double fval = NAN;
    std::vector<double> gradient(model.numParameters(), NAN);

    model.evaluate(parameters, 1.0, fval, gradient, nullptr, nullptr);
    CHECK_EQUAL(3.0, fval);
    CHECK_EQUAL(1.0, gradient[0]);
    CHECK_EQUAL(0.0, gradient[1]);

    std::vector<double> dataset {2.0, 3.0};
    model.evaluate(parameters, dataset, fval, gradient, nullptr, nullptr);
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

    gradFun.evaluate(parameters, fval, gradient);
    CHECK_EQUAL(9.0, fval);
    CHECK_EQUAL(5.0, gradient[0]);
    CHECK_EQUAL(0.0, gradient[1]);
}


#ifdef PARPE_ENABLE_IPOPT
TEST(optimizationProblem, linearModelToGradientFunOptimization) {
    // create optimization problem for the linear model
    // does not do anything meaningful yet
    std::vector<double> dataset {2.0, 3.0};
    auto model = std::unique_ptr<parpe::SummedGradientFunction<double>>(new SummedGradientFunctionLinearModelTest());
    auto gradFun = std::unique_ptr<parpe::GradientFunction>(
                new parpe::SummedGradientFunctionGradientFunctionAdapter<double>(
                    std::move(model), dataset));
    parpe::OptimizationProblemImpl problem {std::move(gradFun), std::make_unique<parpe::Logger>()};

    std::vector<double> parametersMin {-100, -100};
    std::vector<double> parametersMax {100, 100};

    problem.setParametersMin(parametersMin);
    problem.setParametersMax(parametersMax);

    parpe::OptimizerIpOpt opt;
    opt.optimize(&problem);
}
#endif
