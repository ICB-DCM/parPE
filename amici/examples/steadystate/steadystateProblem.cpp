#include "steadystateProblem.h"
#include "wrapfunctions.h"
#include <cstring>

SteadystateProblem::SteadystateProblem()
{
    setupUserData();
    setupExpData();

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

    if(objFunGrad) {
        udata->sensi = AMI_SENSI_ORDER_FIRST;
        udata->sensi_meth = AMI_SENSI_FSA;
    } else {
        udata->sensi = AMI_SENSI_ORDER_NONE;
        udata->sensi_meth = AMI_SENSI_NONE;
    }

    ReturnData *rdata = getSimulationResults(udata, edata);
    int status = (int) *rdata->status;

    *objFunVal = - *rdata->llh;

    if(objFunGrad)
        for(int i = 0; i < udata->np; ++i)
            objFunGrad[i] = - rdata->sllh[i];

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
    printf("Optimal parameters:\n\t");
    printArray(optimalParameters, udata->np);
    printf("\n");
    printf("Minimal cost: %f\n", optimalCost);
}

SteadystateProblem::~SteadystateProblem(){
    delete[] initialParameters;
    delete[] parametersMin;
    delete[] parametersMax;
    delete udata;
    freeExpData(edata);

    delete optimizationOptions;
}

void SteadystateProblem::setupUserData()
{
    udata = new UserData(getUserData());

    udata->nt = 1;
    udata->ts = new double[udata->nt];
    udata->ts[0] = 100;

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
    udata->k[0] = 0.1;
    udata->k[1] = 0.4;
    udata->k[2] = 0.7;
    udata->k[3] = 1;

    udata->sensi = AMI_SENSI_ORDER_FIRST;
    udata->sensi_meth = AMI_SENSI_FSA;
}

void SteadystateProblem::setupExpData()
{
    edata = new ExpData();
    edata->am_my = new double[udata->nytrue * udata->nt];
    fillArray(edata->am_my, udata->nytrue * udata->nt, 1);
    edata->am_ysigma = new double[udata->nytrue * udata->nt];
    fillArray(edata->am_ysigma, udata->nytrue * udata->nt, 1);
}

