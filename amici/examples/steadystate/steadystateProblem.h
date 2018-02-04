#ifndef STEADYSTATEPROBLEM_H
#define STEADYSTATEPROBLEM_H

#include "include/amici_interface_cpp.h"
#include "optimizationProblem.h"
#include <H5Cpp.h>
#include <amici.h>


/**
 * @brief Cost function for the AMICI steady-state example
 */
class ExampleSteadystateGradientFunction : public parpe::GradientFunction {
public:
    ExampleSteadystateGradientFunction(hid_t fileId);
    parpe::FunctionEvaluationStatus evaluate(
            const double* const parameters,
            double &fval,
            double* gradient) const override;

    int numParameters() const override;
    void setupUserData(int conditionIdx);
    void setupExpData(int conditionIdx);

private:
    void requireSensitivities(bool sensitivitiesRequired) const;
    void readFixedParameters(int conditionIdx) const;
    void readMeasurement(int conditionIdx) const;

    hid_t fileId = -1;

    std::unique_ptr<amici::ExpData> edata;
    std::unique_ptr<amici::Model> model;
    std::unique_ptr<amici::Solver> solver;
};


/**
 * @brief Optimization problem for the AMICI steady-state example
 */

class ExampleSteadystateProblem : public parpe::OptimizationProblem {
  public:
    ExampleSteadystateProblem(std::string const& dataFileName);

    void fillInitialParameters(double *buffer) const override;
    void fillParametersMin(double *buffer) const override;
    void fillParametersMax(double *buffer) const override;

private:
    H5::H5File file;
};

#endif // STEADYSTATEPROBLEM_H
