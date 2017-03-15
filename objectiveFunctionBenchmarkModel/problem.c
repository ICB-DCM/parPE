#include "problem.h"
#include "include/rdata.h"
#include "objectiveFunction.h"
#include <math.h>
#include <string.h>


int objectiveFunctionWrapper(void *problem, const double *parameters, double *result);

int objectiveFunctionGradientWrapper(void *problem, const double *parameters, double *objFunVal, double *objFunGrad);

int intermediateFunction(void *_problem,
                          int alg_mod,
                          int iter_count,
                          double obj_value,
                          double inf_pr, double inf_du,
                          double mu,
                          double d_norm,
                          double regularization_size,
                          double alpha_du, double alpha_pr,
                          int ls_trials);

OptimizationProblem *getBenchmarkOptimizationProblem(Datapath path) {
    OptimizationProblem *problem = malloc(sizeof(*problem));
    memset(problem, 0, sizeof(*problem));

    problem->numOptimizationParameters = getLenTheta();
    problem->maxOptimizerIterations = getMaxIter();

    problem->objectiveFunction = objectiveFunctionWrapper;
    problem->objectiveFunctionGradient = objectiveFunctionGradientWrapper;
    problem->intermediateFunction = intermediateFunction;

    MyUserData *myUserData = malloc(sizeof *myUserData);
    myUserData->nTheta   = problem->numOptimizationParameters;
    myUserData->gradient = malloc(sizeof(double) * myUserData->nTheta);
    myUserData->theta    = malloc(sizeof(double) * myUserData->nTheta);
    myUserData->datapath = path;
    myUserData->scaling  = AMI_SCALING_LOG10;
    problem->userData = myUserData;

    problem->initialParameters = malloc(sizeof(double) * problem->numOptimizationParameters);
    getRandomInitialThetaFromFile(path, problem->initialParameters, (AMI_parameter_scaling) myUserData->scaling);

    problem->parametersMin = malloc(sizeof(double) * problem->numOptimizationParameters);
    getThetaLowerBounds(path, problem->parametersMin, myUserData->scaling);
    problem->parametersMax = malloc(sizeof(double) * problem->numOptimizationParameters);
    getThetaUpperBounds(path, problem->parametersMax, myUserData->scaling);


    problem->logOptimizerFinished = saveLocalOptimizerResults;
    problem->logObjectiveFunctionEvaluation = logLocalOptimizerObjectiveFunctionEvaluation;
    problem->logObjectiveFunctionGradientEvaluation = logLocalOptimizerObjectiveFunctionGradientEvaluation;

    return problem;
}

void freeBenchmarkProblem(OptimizationProblem *problem)
{
    MyUserData *myUserData = (MyUserData *) problem->userData;
    free(myUserData->gradient);
    free(myUserData->theta);
    free(myUserData);

    free(problem->initialParameters);
    free(problem->parametersMax);
    free(problem->parametersMin);

    free(problem);
}


int objectiveFunctionWrapper(void *_problem, const double *parameters, double *result) {
//    *result = 1;
//    return 1;

    OptimizationProblem *problem = (OptimizationProblem *) _problem;
    MyUserData *myUserData = (MyUserData *) problem->userData;

    int status = evaluateObjectiveFunction(parameters, problem->numOptimizationParameters,
                                     myUserData->datapath, result, NULL, myUserData->scaling);

    // TODO is this still required?
    myUserData->objectiveFunctionValue = *result;
    for(int i = 0; i < problem->numOptimizationParameters; ++i)
        myUserData->theta[i] = parameters[i];
    fillArray(myUserData->gradient, problem->numOptimizationParameters, NAN);

    return status;
}


int objectiveFunctionGradientWrapper(void *_problem, const double *parameters, double *objFunVal, double *objFunGrad) {
    OptimizationProblem *problem = (OptimizationProblem *) _problem;
    MyUserData *myUserData = (MyUserData *) problem->userData;

//    for(int i = 0; i < problem->numOptimizationParameters; ++i)
//        objFunGrad = 0;
//    *objFunVal = 1;
//    return 1;

    int status = evaluateObjectiveFunction(parameters, problem->numOptimizationParameters,
                                           myUserData->datapath, objFunVal,
                                           objFunGrad, myUserData->scaling);

    myUserData->objectiveFunctionValue = *objFunVal;
    for(int i = 0; i < problem->numOptimizationParameters; ++i) {
        myUserData->theta[i] = parameters[i];
        myUserData->gradient[i] = objFunGrad[i];
    }


    return status;

}


int intermediateFunction(void *_problem,
                          int alg_mod,
                          int iter_count,
                          double obj_value,
                          double inf_pr, double inf_du,
                          double mu,
                          double d_norm,
                          double regularization_size,
                          double alpha_du, double alpha_pr,
                          int ls_trials) {

    OptimizationProblem *problem = (OptimizationProblem *) _problem;
    MyUserData *data = (MyUserData *) problem->userData;
    data->datapath.idxLocalOptimizationIteration = iter_count;

    char strBuf[50];
    sprintDatapath(strBuf, data->datapath);
    //    logmessage(LOGLVL_INFO, "%s: %d %d %e %e %e %e %e %e %e %e %d", strBuf,
    //               alg_mod, iter_count, obj_value, inf_pr, inf_du,
    //               mu, d_norm, regularization_size, alpha_du, alpha_pr, ls_trials);

    logLocalOptimizerIteration(data->datapath, iter_count, data->theta, obj_value, data->gradient, 0, data->nTheta,
                               alg_mod, inf_pr, inf_du, mu, d_norm, regularization_size, alpha_du, alpha_pr, ls_trials);
    return true;
}
