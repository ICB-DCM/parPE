#include "SteadyStateMultiConditionProblem.h"


SteadyStateMultiConditionDataProvider::SteadyStateMultiConditionDataProvider(const char *hdf5Filename)
    : MultiConditionDataProvider(hdf5Filename){
    udata = new UserData(getModelDims());
    setupUserData(udata);
    modelDims.nt = 20;
}

int SteadyStateMultiConditionDataProvider::updateFixedSimulationParameters(int conditionIdx, UserData *udata) {
    hdf5Read2DDoubleHyperslab(fileId, "/data/k", udata->nk, 1, 0, conditionIdx, udata->k);
    return 0;
}

ExpData *SteadyStateMultiConditionDataProvider::getExperimentalDataForCondition(int conditionIdx) {
    ExpData *edata = new ExpData();

    edata->am_my = new double[modelDims.nytrue * modelDims.nt];
    hdf5Read3DDoubleHyperslab(fileId, "/data/ymeasured", 1, modelDims.ny, modelDims.nt, conditionIdx, 0, 0, edata->am_my);

    double ysigma = AMI_HDF5_getDoubleScalarAttribute(fileId, "data", "sigmay");
    edata->am_ysigma = new double[modelDims.nytrue * modelDims.nt];
    fillArray(edata->am_ysigma, modelDims.nytrue * modelDims.nt, ysigma);

    return edata;
}

void SteadyStateMultiConditionDataProvider::setupUserData(UserData *udata) {

    udata->nt = 20;

    hsize_t length;
    AMI_HDF5_getDoubleArrayAttribute(fileId, "data", "t", &udata->ts, &length);
    assert(length == (unsigned) udata->nt);

    udata->idlist = new double[udata->nx];
    fillArray(udata->idlist, udata->nx, 1);
    udata->qpositivex = new double[udata->nx];
    fillArray(udata->qpositivex, udata->nx, 1);

    // calculate sensitivities for all parameters
    udata->plist = new int[udata->np];
    udata->nplist = udata->np;
    for(int i = 0; i < udata->np; ++i) udata->plist[i] = i;

    udata->p = new double[udata->np];

    // set model constants
    udata->k = new double[udata->nk];
    updateFixedSimulationParameters(0, udata);

    udata->maxsteps = 1e5;

    udata->sensi = AMI_SENSI_ORDER_FIRST;
    udata->sensi_meth = AMI_SENSI_FSA;

}

SteadyStateMultiConditionProblem::SteadyStateMultiConditionProblem(SteadyStateMultiConditionDataProvider *dp) : MultiConditionProblem(dp) {

    udata = new UserData(getUserData());
    dp->setupUserData(udata);
    numOptimizationParameters = udata->np;

    initialParameters = new double [numOptimizationParameters];
    fillArray(initialParameters, udata->np, 0);

    parametersMin = new double [numOptimizationParameters];
    fillArray(parametersMin, udata->np, -5);

    parametersMax = new double [numOptimizationParameters];
    fillArray(parametersMax, udata->np, 5);


    optimizationOptions = new OptimizationOptions();
    optimizationOptions->optimizer = OPTIMIZER_IPOPT;
    optimizationOptions->printToStdout = true;
    optimizationOptions->maxOptimizerIterations = 30;

}

void SteadyStateMultiConditionProblem::setSensitivityOptions(bool sensiRequired) {
    // sensitivities requested?
    if(sensiRequired) {
        udata->sensi = AMI_SENSI_ORDER_FIRST;
        udata->sensi_meth = AMI_SENSI_FSA;
    } else {
        udata->sensi = AMI_SENSI_ORDER_NONE;
        udata->sensi_meth = AMI_SENSI_NONE;
    }

}
