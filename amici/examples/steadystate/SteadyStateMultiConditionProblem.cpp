#include "SteadyStateMultiConditionProblem.h"
#include "optimizationOptions.h"
#include <amici_model.h>
#include <cstdio>
#include <misc.h>
#include <multiConditionProblemResultWriter.h>
Model *getModel();

SteadyStateMultiConditionDataProvider::SteadyStateMultiConditionDataProvider(
    Model *model, std::string hdf5Filename)
    : MultiConditionDataProvider(model, hdf5Filename) {

    // hdf5MeasurementPath = "/data/ytrue";

    hdf5MeasurementPath = "/data/ymeasured";
    hdf5MeasurementSigmaPath = "/data/sigmay";

    udata = std::unique_ptr<UserData>(getUserData());
}

int SteadyStateMultiConditionDataProvider::updateFixedSimulationParameters(int conditionIdx, UserData &udata) const {
    parpe::hdf5Read2DDoubleHyperslab(fileId, "/data/k", model->nk, 1, 0, conditionIdx,
                              udata.k);
    return 0;
}


void SteadyStateMultiConditionDataProvider::setupUserData(
    UserData *udata) const {

    hsize_t length;
    AMI_HDF5_getDoubleArrayAttribute(fileId, "data", "t", &udata->ts, &length);
    udata->nt = length;

    // calculate sensitivities for all parameters
    udata->requireSensitivitiesForAllParameters();

    // set model constants
    updateFixedSimulationParameters(0, *udata);

    udata->pscale = AMICI_SCALING_LOG10;
    udata->sensi = AMICI_SENSI_ORDER_FIRST;
    udata->sensi_meth = AMICI_SENSI_FSA;
    udata->maxsteps = 1e4;
    udata->newton_maxlinsteps = 100;
    udata->newton_maxsteps = 40;
}

std::unique_ptr<UserData> SteadyStateMultiConditionDataProvider::getUserData() const {
    auto udata = std::unique_ptr<amici::UserData>(model->getNewUserData());
    setupUserData(udata.get());
    return udata;
}


SteadyStateMultiConditionProblem::SteadyStateMultiConditionProblem(
    SteadyStateMultiConditionDataProvider *dp, parpe::LoadBalancerMaster *loadBalancer)
    : MultiConditionProblem(dp, loadBalancer) {

    std::unique_ptr<parpe::OptimizationOptions> options(parpe::OptimizationOptions::fromHDF5(
                                                             dataProvider->getHdf5FileId()));

    optimizationOptions = *options.get();
    std::fill(initialParameters_.begin(), initialParameters_.end(), 0);
    std::fill(parametersMin_.begin(), parametersMin_.end(), -5);
    std::fill(parametersMax_.begin(), parametersMax_.end(), 5);
}
