#include "steadyStateMultiConditionDataprovider.h"
#include "optimizationOptions.h"
#include <amici_model.h>
#include <cstdio>
#include <misc.h>
#include <multiConditionProblemResultWriter.h>

SteadyStateMultiConditionDataProvider::SteadyStateMultiConditionDataProvider(Model *model, std::string hdf5Filename, std::string rootPath)
    : MultiConditionDataProvider(model, hdf5Filename, rootPath) {

    udata = std::unique_ptr<UserData>(getUserData());
}

int SteadyStateMultiConditionDataProvider::getNumConditionSpecificParametersPerSimulation() const {
    return 0;
}


void SteadyStateMultiConditionDataProvider::setupUserData(
    UserData *udata) const {

    hsize_t length;
    auto timePath = rootPath + "/parameters";
    AMI_HDF5_getDoubleArrayAttribute(fileId, timePath.c_str(), "t", &udata->ts, &length);
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
