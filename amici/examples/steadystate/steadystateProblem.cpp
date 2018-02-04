#include "steadystateProblem.h"
#include "hdf5Misc.h"
#include "optimizationOptions.h"
#include "wrapfunctions.h"
#include <amici_hdf5.h>
#include <amici_model.h>
#include <cassert>
#include <cstring>
#include <iostream>
#include <misc.h>

ExampleSteadystateProblem::ExampleSteadystateProblem(const std::string &dataFileName)
    : file(H5::H5File(dataFileName, H5F_ACC_RDONLY))
{
    auto optimizationOptions = getOptimizationOptions();
    optimizationOptions.optimizer = parpe::OPTIMIZER_IPOPT;
    optimizationOptions.printToStdout = true;
    optimizationOptions.maxOptimizerIterations = 100;
    setOptimizationOptions(optimizationOptions);

    costFun = std::make_unique<ExampleSteadystateGradientFunction>(file.getId());
}

void ExampleSteadystateProblem::fillInitialParameters(double *buffer) const
{
    std::fill(buffer, buffer + costFun->numParameters(), 0.0);

}

void ExampleSteadystateProblem::fillParametersMin(double *buffer) const
{
    std::fill(buffer, buffer + costFun->numParameters(), -5.0);

}

void ExampleSteadystateProblem::fillParametersMax(double *buffer) const
{
    std::fill(buffer, buffer + costFun->numParameters(), 5.0);

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
    hsize_t length;
    double *buf;
    amici::AMI_HDF5_getDoubleArrayAttribute(fileId, "parameters", "t", &buf, &length);
    model->setTimepoints(std::vector<double>(buf, buf+length));
    delete[] buf;

    // set model constants
    readFixedParameters(conditionIdx);

    model->setParameterScale(amici::AMICI_SCALING_LOG10);
    requireSensitivities(true);
    model->requireSensitivitiesForAllParameters();
    solver->setMaxSteps(10000);
}

void ExampleSteadystateGradientFunction::setupExpData(int conditionIdx) {
    edata.reset(new amici::ExpData(model.get()));
    readMeasurement(conditionIdx);
}


void ExampleSteadystateGradientFunction::readFixedParameters(int conditionIdx) const {
    std::vector<double> k(model->nk());
    parpe::hdf5Read2DDoubleHyperslab(fileId, "/fixedParameters/k", k.size(), 1, 0, conditionIdx, k.data());
    model->setFixedParameters(k);
}

void ExampleSteadystateGradientFunction::readMeasurement(int conditionIdx) const {
    parpe::hdf5Read3DDoubleHyperslab(fileId, "/measurements/y", 1, model->ny,
                              model->nt(), conditionIdx, 0, 0, edata->my.data());


    parpe::hdf5Read3DDoubleHyperslab(fileId, "/measurements/ysigma", 1, edata->nytrue,
                                     edata->nt, conditionIdx, 0, 0, edata->sigmay.data());
}

ExampleSteadystateGradientFunction::ExampleSteadystateGradientFunction(hid_t fileId)
    : fileId(fileId), model(getModel()), solver(model->getSolver())
{
    setupUserData(0);
    setupExpData(0);
}

parpe::FunctionEvaluationStatus ExampleSteadystateGradientFunction::evaluate(const double * const parameters, double &fval, double *gradient) const
{

    model->setParameters(std::vector<double>(parameters, parameters + numParameters()));

    //    printArray(parameters, udata->np);printf("\n");

    requireSensitivities(gradient);

    std::unique_ptr<amici::ReturnData> rdata {amici::getSimulationResults(*model, edata.get(), *solver)};
    int status = (int)*rdata->status;

    fval = -*rdata->llh;

    if (gradient)
        for (int i = 0; i < model->np(); ++i)
            gradient[i] = -rdata->sllh[i];

    return status == 0 ? parpe::functionEvaluationSuccess : parpe::functionEvaluationFailure;
}

int ExampleSteadystateGradientFunction::numParameters() const
{
    return model->np();
}
