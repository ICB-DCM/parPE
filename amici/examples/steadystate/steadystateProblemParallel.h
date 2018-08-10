#ifndef STEADYSTATEPROBLEM_PARALLEL_H
#define STEADYSTATEPROBLEM_PARALLEL_H

#include "steadyStateMultiConditionDataprovider.h"
#include <loadBalancerMaster.h>
#include <loadBalancerWorker.h>
#include <memory>

#include <gsl/gsl-lite.hpp>

/**
 * @brief The ExampleSteadystateGradientFunctionParallel class evaluates an ODE-constrained objective function in paralell.
 */

class ExampleSteadystateGradientFunctionParallel : public parpe::GradientFunction {
public:
    ExampleSteadystateGradientFunctionParallel(parpe::LoadBalancerMaster *loadBalancer, const std::string &dataFileName);

    parpe::FunctionEvaluationStatus evaluate(
            gsl::span<double const> parameters,
            double &fval,
            gsl::span<double> gradient) const override;

    int numParameters() const override;
    void setupUserData(int conditionIdx);
    void setupExpData(int conditionIdx);

    void messageHandler(std::vector<char> &buffer, int jobId);


private:

    int evaluateParallel(gsl::span<const double> parameters, double &objFunVal,
                         gsl::span<double> objFunGrad) const;

    int evaluateSerial(gsl::span<double const> parameters, double &objFunVal,
                       gsl::span<double> objFunGrad) const;


    void requireSensitivities(bool sensitivitiesRequired) const;
    void readFixedParameters(int conditionIdx) const;
    void readMeasurement(int conditionIdx) const;

    parpe::LoadBalancerMaster *loadBalancer = nullptr; // non-owning

    std::unique_ptr<amici::ExpData> edata;
    std::unique_ptr<amici::Model> model;
    std::unique_ptr<amici::Solver> solver;
    int numConditions;

    std::unique_ptr<SteadyStateMultiConditionDataProvider> dataProvider;
};

#endif // STEADYSTATEPROBLEM_PARALLEL_H
