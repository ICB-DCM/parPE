#include "steadystateProblem.h"
#include "wrapfunctions.h"
#include <cstring>
#include "hdf5Misc.h"
#include "include/ami_hdf5.h"
#include <cassert>

SteadystateProblem::SteadystateProblem()
{
    fileId = H5Fopen("/home/dweindl/src/parPE/amici/examples/steadystate/data.h5", H5F_ACC_RDONLY, H5P_DEFAULT);

    setupUserData(0);
    setupExpData(0);

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
    optimizationOptions->maxOptimizerIterations = 100;
}

int SteadystateProblem::evaluateObjectiveFunction(const double *parameters, double *objFunVal, double *objFunGrad)
{
    memcpy(udata->p, parameters, udata->np * sizeof(double));

//    printArray(parameters, udata->np);printf("\n");

    requireSensitivities(objFunGrad);

    ReturnData *rdata = getSimulationResults(udata, edata);
    int status = (int) *rdata->status;

    *objFunVal = -*rdata->llh;

    if(objFunGrad)
        for(int i = 0; i < udata->np; ++i)
            objFunGrad[i] = -rdata->sllh[i];

    delete rdata;

    return status;
}

int SteadystateProblem::intermediateFunction(int alg_mod, int iter_count, double obj_value, double inf_pr, double inf_du, double mu, double d_norm, double regularization_size, double alpha_du, double alpha_pr, int ls_trials)
{
    return 0;
}

void SteadystateProblem::logObjectiveFunctionEvaluation(const double *parameters, double objectiveFunctionValue, const double *objectiveFunctionGradient, int numFunctionCalls, double timeElapsed)
{
}

void SteadystateProblem::logOptimizerFinished(double optimalCost, const double *optimalParameters, double masterTime, int exitStatus)
{
    printf("Minimal cost: %f\n", optimalCost);
    printf("Optimal parameters  : ");
    printArray(optimalParameters, udata->np);
    printf("\n");
    printf("True parameters were: ");

    hsize_t length;
    double *ptrue;
    AMI_HDF5_getDoubleArrayAttribute(fileId, "/data/", "ptrue", &ptrue, &length);
    printArray(ptrue, length);
    delete ptrue;
    printf("\n");

    printf("Wall time (min): %f\n", masterTime / 60);
}

SteadystateProblem::~SteadystateProblem(){
    H5Fclose(fileId);
    delete[] initialParameters;
    delete[] parametersMin;
    delete[] parametersMax;
    delete udata;
    freeExpData(edata);

    delete optimizationOptions;
}

void SteadystateProblem::requireSensitivities(bool sensitivitiesRequired)
{
    if(sensitivitiesRequired) {
        udata->sensi = AMI_SENSI_ORDER_FIRST;
        udata->sensi_meth = AMI_SENSI_FSA;
    } else {
        udata->sensi = AMI_SENSI_ORDER_NONE;
        udata->sensi_meth = AMI_SENSI_NONE;
    }
}

void SteadystateProblem::setupUserData(int conditionIdx)
{
    udata = new UserData(getUserData());

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
    readFixedParameters(conditionIdx);

    udata->maxsteps = 1e5;

    requireSensitivities(true);
}

void SteadystateProblem::setupExpData(int conditionIdx)
{
    edata = new ExpData();

    edata->am_my = new double[udata->nytrue * udata->nt];
    readMeasurement(conditionIdx);

    double ysigma = AMI_HDF5_getDoubleScalarAttribute(fileId, "data", "sigmay");
    edata->am_ysigma = new double[udata->nytrue * udata->nt];
    fillArray(edata->am_ysigma, udata->nytrue * udata->nt, ysigma);
}

void SteadystateProblem::readFixedParameters(int conditionIdx)
{
    hdf5Read2DDoubleHyperslab(fileId, "/data/k", udata->nk, 1, 0, conditionIdx, udata->k);
}

void SteadystateProblem::readMeasurement(int conditionIdx)
{
    hdf5Read3DDoubleHyperslab(fileId, "/data/ymeasured", 1, udata->ny, udata->nt, conditionIdx, 0, 0, edata->am_my);
}

