#include "SteadyStateMultiConditionProblem.h"
#include "optimizationOptions.h"
#include <amici_model.h>
#include <cstdio>
#include <misc.h>
Model *getModel();

SteadyStateMultiConditionDataProvider::SteadyStateMultiConditionDataProvider(
    Model *model, std::string hdf5Filename)
    : MultiConditionDataProvider(model, hdf5Filename) {

    // hdf5MeasurementPath = "/data/ytrue";

    hdf5MeasurementPath = "/data/ymeasured";
    hdf5MeasurementSigmaPath = "/data/sigmay";

    udata = getUserData();
}

int SteadyStateMultiConditionDataProvider::updateFixedSimulationParameters(
    int conditionIdx, UserData *udata) const {
    parPE::hdf5Read2DDoubleHyperslab(fileId, "/data/k", model->nk, 1, 0, conditionIdx,
                              udata->k);
    return 0;
}


void SteadyStateMultiConditionDataProvider::setupUserData(
    UserData *udata) const {

    hsize_t length;
    AMI_HDF5_getDoubleArrayAttribute(fileId, "data", "t", &udata->ts, &length);
    udata->nt = length;
    udata->qpositivex = new double[model->nx];
    fillArray(udata->qpositivex, model->nx, 1);

    // calculate sensitivities for all parameters
    udata->requireSensitivitiesForAllParameters();
    udata->p = new double[model->np];

    // set model constants
    udata->k = new double[model->nk];
    updateFixedSimulationParameters(0, udata);

    udata->pscale = AMICI_SCALING_LOG10;
    udata->sensi = AMICI_SENSI_ORDER_FIRST;
    udata->sensi_meth = AMICI_SENSI_FSA;

    udata->maxsteps = 1e4;
}

UserData *SteadyStateMultiConditionDataProvider::getUserData() const {
    UserData *udata = model->getNewUserData();
    setupUserData(udata);
    return udata;
}

SteadyStateMultiConditionDataProvider::
    ~SteadyStateMultiConditionDataProvider() {
    delete udata;
}

SteadyStateMultiConditionProblem::SteadyStateMultiConditionProblem(
    SteadyStateMultiConditionDataProvider *dp, parPE::LoadBalancerMaster *loadBalancer)
    : MultiConditionProblem(dp, loadBalancer) {

    std::fill(initialParameters_.begin(), initialParameters_.end(), 0);
    std::fill(parametersMin_.begin(), parametersMin_.end(), -5);
    std::fill(parametersMax_.begin(), parametersMax_.end(), 5);


    optimizationOptions = new parPE::OptimizationOptions();
    optimizationOptions->optimizer = parPE::OPTIMIZER_IPOPT;
    optimizationOptions->printToStdout = true;
    optimizationOptions->maxOptimizerIterations = 30;
}

void SteadyStateMultiConditionProblem::setSensitivityOptions(
    bool sensiRequired) {
    // sensitivities requested?
    if (sensiRequired) {
        udata->sensi = AMICI_SENSI_ORDER_FIRST;
        udata->sensi_meth = AMICI_SENSI_FSA;
    } else {
        udata->sensi = AMICI_SENSI_ORDER_NONE;
        udata->sensi_meth = AMICI_SENSI_NONE;
    }
}

SteadyStateMultiConditionProblem::~SteadyStateMultiConditionProblem() {
    delete optimizationOptions;
}
