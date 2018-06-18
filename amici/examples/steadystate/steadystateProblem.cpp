#include "steadystateProblem.h"
#include "hdf5Misc.h"
#include "optimizationOptions.h"
#include "wrapfunctions.h"
#include <cassert>
#include <cstring>
#include <iostream>
#include <misc.h>
#include <amici/hdf5.h>
#include <hdf5Misc.h>

ExampleSteadystateProblem::ExampleSteadystateProblem(const std::string &dataFileName)
{
    auto lock = parpe::hdf5MutexGetLock();
    file.openFile(dataFileName, H5F_ACC_RDONLY);

    auto optimizationOptions = getOptimizationOptions();
    optimizationOptions.optimizer = parpe::optimizerName::OPTIMIZER_IPOPT;
    optimizationOptions.printToStdout = true;
    optimizationOptions.maxOptimizerIterations = 100;
    setOptimizationOptions(optimizationOptions);

    costFun = std::make_unique<ExampleSteadystateGradientFunction>(file.getId());
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
        solver->setSensitivityOrder(amici::AMICI_SENSI_ORDER_FIRST);
        solver->setSensitivityMethod(amici::AMICI_SENSI_FSA);
    } else {
        solver->setSensitivityOrder(amici::AMICI_SENSI_ORDER_NONE);
        solver->setSensitivityMethod(amici::AMICI_SENSI_NONE);
    }
}

void ExampleSteadystateGradientFunction::setupUserData(int conditionIdx) {
    hsize_t m = 0, n = 0;
    model->setTimepoints(amici::hdf5::getDoubleDataset2D(fileId, "/parameters/t", m, n));

    // set model constants
    readFixedParameters(conditionIdx);

    model->setParameterScale(amici::AMICI_SCALING_LOG10);
    requireSensitivities(true);
    model->requireSensitivitiesForAllParameters();
    solver->setMaxSteps(10000);
}

void ExampleSteadystateGradientFunction::setupExpData(int conditionIdx) {
    edata.reset(new amici::ExpData(*model));
    readMeasurement(conditionIdx);
}


void ExampleSteadystateGradientFunction::readFixedParameters(int conditionIdx) const {
    std::vector<double> k(model->nk());
    parpe::hdf5Read2DDoubleHyperslab(fileId, "/fixedParameters/k", k.size(), 1, 0, conditionIdx, k.data());
    model->setFixedParameters(k);
}

void ExampleSteadystateGradientFunction::readMeasurement(int conditionIdx) const {
    parpe::hdf5Read3DDoubleHyperslab(fileId, "/measurements/y",
                                     1, edata->nt, edata->nytrue,
                                     conditionIdx, 0, 0,
                                     edata->my.data());

    parpe::hdf5Read3DDoubleHyperslab(fileId, "/measurements/ysigma",
                                     1, edata->nt, edata->nytrue,
                                     conditionIdx, 0, 0,
                                     edata->sigmay.data());
}

ExampleSteadystateGradientFunction::ExampleSteadystateGradientFunction(hid_t fileId)
    : fileId(fileId), model(getModel()), solver(model->getSolver())
{
    setupUserData(0);
    setupExpData(0);
}

parpe::FunctionEvaluationStatus ExampleSteadystateGradientFunction::evaluate(gsl::span<const double> parameters, double &fval, gsl::span<double> gradient) const
{

    model->setParameters(std::vector<double>(parameters.begin(), parameters.end()));

    //    printArray(parameters, udata->np);printf("\n");

    requireSensitivities(gradient.size());

    auto rdata = amici::runAmiciSimulation(*solver, edata.get(), *model);

    fval = -rdata->llh;

    if (gradient.size())
        for (int i = 0; i < model->np(); ++i)
            gradient[i] = -rdata->sllh[i];

    return rdata->status == 0 ? parpe::functionEvaluationSuccess : parpe::functionEvaluationFailure;
}

int ExampleSteadystateGradientFunction::numParameters() const
{
    return model->np();
}
