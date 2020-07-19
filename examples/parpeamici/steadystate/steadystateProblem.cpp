#include "steadystateProblem.h"

#include <parpecommon/hdf5Misc.h>
#include <parpeoptimization/optimizationOptions.h>
#include <parpecommon/misc.h>

#include "wrapfunctions.h"

#include <cassert>
#include <cstring>
#include <iostream>

#include <amici/hdf5.h>

ExampleSteadystateProblem::ExampleSteadystateProblem(const std::string &dataFileName)
{
    auto lock = parpe::hdf5MutexGetLock();
    file.openFile(dataFileName, H5F_ACC_RDONLY);

    auto optimizationOptions = OptimizationProblem::getOptimizationOptions();
    optimizationOptions.optimizer = parpe::optimizerName::OPTIMIZER_IPOPT;
    optimizationOptions.printToStdout = true;
    optimizationOptions.maxOptimizerIterations = 100;
    OptimizationProblem::setOptimizationOptions(optimizationOptions);

    cost_fun_ = std::make_unique<ExampleSteadystateGradientFunction>(file.getId());
}

void ExampleSteadystateProblem::fillInitialParameters(gsl::span<double> buffer) const
{
    std::fill(buffer.begin(), buffer.end(), 0.0);

}

void ExampleSteadystateProblem::fillParametersMin(gsl::span<double> buffer) const
{
    std::fill(buffer.begin(), buffer.end(), -5.0);

}

void ExampleSteadystateProblem::fillParametersMax(gsl::span<double> buffer) const
{
    std::fill(buffer.begin(), buffer.end(), 5.0);

}

void ExampleSteadystateGradientFunction::requireSensitivities(
    bool sensitivitiesRequired) const {
    if (sensitivitiesRequired) {
        solver->setSensitivityOrder(amici::SensitivityOrder::first);
        solver->setSensitivityMethod(amici::SensitivityMethod::forward);
    } else {
        solver->setSensitivityOrder(amici::SensitivityOrder::none);
        solver->setSensitivityMethod(amici::SensitivityMethod::none);
    }
}

void ExampleSteadystateGradientFunction::setupUserData(int conditionIdx) {
    hsize_t m = 0, n = 0;
    auto lock = parpe::hdf5MutexGetLock();
    model->setTimepoints(amici::hdf5::getDoubleDataset2D(fileId, "/parameters/t", m, n));

    // set model constants
    readFixedParameters(conditionIdx);

    model->setParameterScale(amici::ParameterScaling::log10);
    requireSensitivities(true);
    model->requireSensitivitiesForAllParameters();
    solver->setMaxSteps(10000);
}

void ExampleSteadystateGradientFunction::setupExpData(int conditionIdx) {
    edata = std::make_unique<amici::ExpData>(*model);
    readMeasurement(conditionIdx);
}

std::vector<std::string> ExampleSteadystateGradientFunction::getParameterIds() const
{
    return parpe::hdf5Read1dStringDataset(fileId, "/parameters/parameterNames");
}


void ExampleSteadystateGradientFunction::readFixedParameters(int conditionIdx) const {
    std::vector<double> k(model->nk());
    parpe::hdf5Read2DDoubleHyperslab(fileId, "/fixedParameters/k", k.size(),
                                     1, 0, conditionIdx, k);
    model->setFixedParameters(k);
}

void ExampleSteadystateGradientFunction::readMeasurement(int conditionIdx) const {
    edata->setObservedData(
                parpe::hdf5Get3DDoubleHyperslab(fileId, "/measurements/y",
                                                1, edata->nt(), edata->nytrue(),
                                                conditionIdx, 0, 0));
    edata->setObservedDataStdDev(
                parpe::hdf5Get3DDoubleHyperslab(fileId, "/measurements/ysigma",
                                                1, edata->nt(), edata->nytrue(),
                                                conditionIdx, 0, 0));
}

ExampleSteadystateGradientFunction::ExampleSteadystateGradientFunction(hid_t fileId)
    : fileId(fileId), model(amici::generic_model::getModel()), solver(model->getSolver())
{
    setupUserData(0);
    setupExpData(0);
}

parpe::FunctionEvaluationStatus ExampleSteadystateGradientFunction::evaluate(
        gsl::span<const double> parameters, double &fval,
        gsl::span<double> gradient, parpe::Logger * /*logger*/,
        double * /*cpuTime*/) const
{

    model->setParameters(std::vector<double>(parameters.begin(), parameters.end()));

    //    printArray(parameters, udata->np);printf("\n");

    requireSensitivities(!gradient.empty());

    auto rdata = amici::runAmiciSimulation(*solver, edata.get(), *model);

    fval = -rdata->llh;

    if (!gradient.empty())
        for (int i = 0; i < model->np(); ++i)
            gradient[i] = -rdata->sllh[i];

    return rdata->status == 0 ? parpe::functionEvaluationSuccess : parpe::functionEvaluationFailure;
}

int ExampleSteadystateGradientFunction::numParameters() const
{
    return model->np();
}
