#include "steadystateProblem.h"
#include "hdf5Misc.h"
#include "optimizationOptions.h"
#include "wrapfunctions.h"
#include <amici_hdf5.h>
#include <amici_model.h>
#include <cassert>
#include <cstring>
#include <iostream>

#define STEADYSTATE_DATA_HDF5 "/home/dweindl/src/parPE/amici/examples/steadystate/data.h5"

ExampleSteadystateProblem::ExampleSteadystateProblem() {
    fileId = H5Fopen(STEADYSTATE_DATA_HDF5, H5F_ACC_RDONLY, H5P_DEFAULT);

    model = getModel();
    setupUserData(0);
    setupExpData(0);

    setNumOptimizationParameters(model->np);
    fillArray(initialParameters_.data(), model->np, 0);
    fillArray(parametersMin_.data(), model->np, -5);
    fillArray(parametersMax_.data(), model->np, 5);

    optimizationOptions = new OptimizationOptions();
    optimizationOptions->optimizer = OPTIMIZER_IPOPT;
    optimizationOptions->printToStdout = true;
    optimizationOptions->maxOptimizerIterations = 100;
}

int ExampleSteadystateProblem::evaluateObjectiveFunction(const double *parameters, double *objFunVal, double *objFunGrad) {

    udata->setParameters(parameters);

    //    printArray(parameters, udata->np);printf("\n");

    requireSensitivities(objFunGrad);

    ReturnData *rdata = getSimulationResults(model, udata, edata);
    int status = (int)*rdata->status;

    *objFunVal = -*rdata->llh;

    if (objFunGrad)
        for (int i = 0; i < model->np; ++i)
            objFunGrad[i] = -rdata->sllh[i];

    delete rdata;

    return status;
}

int ExampleSteadystateProblem::intermediateFunction(
    int alg_mod, int iter_count, double obj_value, double inf_pr, double inf_du,
    double mu, double d_norm, double regularization_size, double alpha_du,
    double alpha_pr, int ls_trials) {
    return 0;
}

void ExampleSteadystateProblem::logObjectiveFunctionEvaluation(
    const double *parameters, double objectiveFunctionValue,
    const double *objectiveFunctionGradient, int numFunctionCalls,
    double timeElapsed) {}

void ExampleSteadystateProblem::logOptimizerFinished(
    double optimalCost, const double *optimalParameters, double masterTime,
    int exitStatus) {
    printf("Minimal cost: %f\n", optimalCost);
    printf("Optimal parameters  : ");
    printArray(optimalParameters, model->np);
    printf("\n");
    printf("True parameters were: ");

    hsize_t length;
    double *ptrue;
    AMI_HDF5_getDoubleArrayAttribute(fileId, "/data/", "ptrue", &ptrue,
                                     &length);
    printArray(ptrue, length);
    delete ptrue;
    printf("\n");

    printf("Wall time (min): %f\n", masterTime / 60);
}

ExampleSteadystateProblem::~ExampleSteadystateProblem() {
    H5Fclose(fileId);
    delete udata;
    delete edata;
    if (model)
        delete model;

    delete optimizationOptions;
}

void ExampleSteadystateProblem::requireSensitivities(
    bool sensitivitiesRequired) {
    if (sensitivitiesRequired) {
        udata->sensi = AMICI_SENSI_ORDER_FIRST;
        udata->sensi_meth = AMICI_SENSI_FSA;
    } else {
        udata->sensi = AMICI_SENSI_ORDER_NONE;
        udata->sensi_meth = AMICI_SENSI_NONE;
    }
}

void ExampleSteadystateProblem::setupUserData(int conditionIdx) {
    udata = model->getNewUserData();

    udata->qpositivex = new double[model->nx];
    fillArray(udata->qpositivex, model->nx, 1);

    hsize_t length;
    AMI_HDF5_getDoubleArrayAttribute(fileId, "data", "t", &udata->ts, &length);
    udata->nt = length;

    // set model constants
    udata->k = new double[model->nk];
    readFixedParameters(conditionIdx);

    udata->pscale = AMICI_SCALING_LOG10;
    requireSensitivities(true);
    udata->requireSensitivitiesForAllParameters();
    udata->maxsteps = 10000;
}

void ExampleSteadystateProblem::setupExpData(int conditionIdx) {
    edata = new ExpData(udata, model);
    readMeasurement(conditionIdx);

    double ysigma = NAN;
    AMI_HDF5_getDoubleScalarAttribute(fileId, "data", "sigmay", &ysigma);
    fillArray(edata->sigmay, model->nytrue * udata->nt, ysigma);
}

void ExampleSteadystateProblem::readFixedParameters(int conditionIdx) {
    hdf5Read2DDoubleHyperslab(fileId, "/data/k", model->nk, 1, 0, conditionIdx,
                              udata->k);
}

void ExampleSteadystateProblem::readMeasurement(int conditionIdx) {
    hdf5Read3DDoubleHyperslab(fileId, "/data/ymeasured", 1, model->ny,
                              udata->nt, conditionIdx, 0, 0, edata->my);
}
