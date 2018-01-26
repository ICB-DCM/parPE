#ifndef STEADYSTATEPROBLEM_PARALLEL_H
#define STEADYSTATEPROBLEM_PARALLEL_H

#include "steadyStateMultiConditionDataprovider.h"
#include <LoadBalancerMaster.h>
#include <LoadBalancerWorker.h>
#include <memory>

/**
 * @brief The ExampleSteadystateGradientFunctionParallel class evaluates an ODE-constrained objective function in paralell.
 */

class ExampleSteadystateGradientFunctionParallel : public parpe::GradientFunction {
public:
    ExampleSteadystateGradientFunctionParallel(parpe::LoadBalancerMaster *loadBalancer, const std::string &dataFileName);

    FunctionEvaluationStatus evaluate(
            const double* const parameters,
            double &fval,
            double* gradient) const override;

    int numParameters() const override;
    void setupUserData(int conditionIdx);
    void setupExpData(int conditionIdx);

    void messageHandler(std::vector<char> &buffer, int jobId);


private:

    int evaluateParallel(const double *parameters, double &objFunVal,
                         double *objFunGrad) const;

    int evaluateSerial(const double *parameters, double &objFunVal,
                       double *objFunGrad) const;


    void requireSensitivities(bool sensitivitiesRequired) const;
    void readFixedParameters(int conditionIdx) const;
    void readMeasurement(int conditionIdx) const;

    parpe::LoadBalancerMaster *loadBalancer = nullptr;

    std::unique_ptr<amici::UserData> udata;
    std::unique_ptr<amici::ExpData> edata;
    std::unique_ptr<amici::Model> model;
    int numConditions;

    std::unique_ptr<SteadyStateMultiConditionDataProvider> dataProvider;
};

#endif // STEADYSTATEPROBLEM_PARALLEL_H
