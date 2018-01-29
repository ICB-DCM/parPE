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
    fillArray(buffer, costFun->numParameters(), 0);

}

void ExampleSteadystateProblem::fillParametersMin(double *buffer) const
{
    fillArray(buffer, costFun->numParameters(), -5);

}

void ExampleSteadystateProblem::fillParametersMax(double *buffer) const
{
    fillArray(buffer, costFun->numParameters(), 5);

}

void ExampleSteadystateGradientFunction::requireSensitivities(
    bool sensitivitiesRequired) const {
    if (sensitivitiesRequired) {
        udata->sensi = AMICI_SENSI_ORDER_FIRST;
        udata->sensi_meth = AMICI_SENSI_FSA;
    } else {
        udata->sensi = AMICI_SENSI_ORDER_NONE;
        udata->sensi_meth = AMICI_SENSI_NONE;
    }
}

void ExampleSteadystateGradientFunction::setupUserData(int conditionIdx) {
    udata.reset(model->getNewUserData());

    hsize_t length;
    AMI_HDF5_getDoubleArrayAttribute(fileId, "parameters", "t", &udata->ts, &length);
    udata->nt = length;

    // set model constants
    readFixedParameters(conditionIdx);

    udata->pscale = AMICI_SCALING_LOG10;
    requireSensitivities(true);
    udata->requireSensitivitiesForAllParameters();
    udata->maxsteps = 10000;
}

void ExampleSteadystateGradientFunction::setupExpData(int conditionIdx) {
    edata.reset(new ExpData(udata.get(), model.get()));
    readMeasurement(conditionIdx);
}


void ExampleSteadystateGradientFunction::readFixedParameters(int conditionIdx) const {
    parpe::hdf5Read2DDoubleHyperslab(fileId, "/fixedParameters/k", model->nk, 1, 0, conditionIdx,
                              udata->k);
}

void ExampleSteadystateGradientFunction::readMeasurement(int conditionIdx) const {
    parpe::hdf5Read3DDoubleHyperslab(fileId, "/measurements/y", 1, model->ny,
                              udata->nt, conditionIdx, 0, 0, edata->my);


    parpe::hdf5Read3DDoubleHyperslab(fileId, "/measurements/ysigma", 1, edata->nytrue,
                                     edata->nt, conditionIdx, 0, 0, edata->sigmay);
}

ExampleSteadystateGradientFunction::ExampleSteadystateGradientFunction(hid_t fileId) : model(getModel()), fileId(fileId)
{
    setupUserData(0);
    setupExpData(0);
}

parpe::FunctionEvaluationStatus ExampleSteadystateGradientFunction::evaluate(const double * const parameters, double &fval, double *gradient) const
{

    udata->setParameters(parameters);

    //    printArray(parameters, udata->np);printf("\n");

    requireSensitivities(gradient);

    std::unique_ptr<amici::ReturnData> rdata {getSimulationResults(model.get(), udata.get(), edata.get())};
    int status = (int)*rdata->status;

    fval = -*rdata->llh;

    if (gradient)
        for (int i = 0; i < model->np; ++i)
            gradient[i] = -rdata->sllh[i];

    return status == 0 ? parpe::functionEvaluationSuccess : parpe::functionEvaluationFailure;
}

int ExampleSteadystateGradientFunction::numParameters() const
{
    return model->np;
}
