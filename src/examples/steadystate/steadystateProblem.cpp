#include "steadystateProblem.h"
#include "wrapfunctions.h"
#include <cstring>

SteadystateProblem::SteadystateProblem()
{
    setupUserData();
    setupExpData();

    numOptimizationParameters = udata->am_np;

    initialParameters = new double [numOptimizationParameters];
    fillArray(initialParameters, udata->am_np, 0);

    parametersMin = new double [numOptimizationParameters];
    fillArray(parametersMin, udata->am_np, -5);

    parametersMax = new double [numOptimizationParameters];
    fillArray(parametersMax, udata->am_np, 5);

    optimizer = OPTIMIZER_IPOPT;

    printToStdout = true;

    maxOptimizerIterations = 100;
}

int SteadystateProblem::evaluateObjectiveFunction(const double *parameters, double *objFunVal, double *objFunGrad)
{
    int status = -1;
    memcpy(udata->am_p, parameters, udata->am_np * sizeof(double));

//    printArray(parameters, udata->am_np);printf("\n");

    if(objFunGrad) {
        udata->am_sensi = AMI_SENSI_ORDER_FIRST;
        udata->am_sensi_meth = AMI_SENSI_FSA;
    } else {
        udata->am_sensi = AMI_SENSI_ORDER_NONE;
        udata->am_sensi_meth = AMI_SENSI_NONE;
    }

    ReturnData *rdata = getSimulationResults(udata, edata, &status);

    *objFunVal = - *rdata->am_llhdata;

    if(objFunGrad)
        for(int i = 0; i < udata->am_np; ++i)
            objFunGrad[i] = - rdata->am_sllhdata[i];

    freeReturnData(rdata);

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
    printArray(optimalParameters, udata->am_np);
    printf("\n");
    printf("Minimal cost: %f\n", optimalCost);
}

SteadystateProblem::~SteadystateProblem(){
    delete[] initialParameters;
    delete[] parametersMin;
    delete[] parametersMax;
    freeUserData(udata);
    freeExpData(edata);
}

void SteadystateProblem::setupUserData()
{
    udata = getDefaultUserData();
    init_modeldims(udata);

    udata->am_atol = 1e-8;
    udata->am_rtol = 1e-8;

    udata->am_nt = 1;
    udata->am_ts = new double[udata->am_nt];
    udata->am_ts[0] = 100;

    udata->am_idlist = new double[udata->am_nx];
    fillArray(udata->am_idlist, udata->am_nx, 1);
    udata->am_qpositivex = new double[udata->am_nx];
    fillArray(udata->am_qpositivex, udata->am_nx, 1);

    udata->am_plist = new int[udata->am_np];
    udata->am_nplist = udata->am_np;
    for(int i = 0; i < udata->am_np; ++i) udata->am_plist[i] = i;

    udata->am_p = new double[udata->am_np];

    udata->am_k = new double[udata->am_nk];
    udata->am_k[0] = 0.1;
    udata->am_k[1] = 0.4;
    udata->am_k[2] = 0.7;
    udata->am_k[3] = 1;

    udata->am_lmm = 1;
    udata->am_iter = 1;
    udata->am_linsol = AMI_KLU;

    udata->am_maxsteps = 1e5;

    udata->am_sensi = AMI_SENSI_ORDER_FIRST;
    udata->am_sensi_meth = AMI_SENSI_FSA;

    processUserData(udata);
}

void SteadystateProblem::setupExpData()
{
    edata = new ExpData();
    edata->am_my = new double[udata->am_nytrue * udata->am_nt];
    fillArray(edata->am_my, udata->am_nytrue * udata->am_nt, 1);
    edata->am_ysigma = new double[udata->am_nytrue * udata->am_nt];
    fillArray(edata->am_ysigma, udata->am_nytrue * udata->am_nt, 1);
}

