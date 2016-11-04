#include "localoptimization.h"

#include <stdbool.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include <signal.h>
#include <alloca.h>

#include "objectivefunction.h"

extern volatile sig_atomic_t caughtTerminationSignal;

typedef struct {
    int nTheta;
    double *theta;
    UserData udata;
    ExpData edata;
    ReturnData rdata;
    double *gradient;
    double objectiveFunctionValue;
    datapath datapath;
} MyUserData;

static IpoptProblem setupIpoptProblem();

/******************************/

static Bool Eval_F(Index n, Number* x, Bool new_x, Number* obj_value, UserDataPtr user_data);

static Bool Eval_Grad_F(Index n, Number* x, Bool new_x, Number* grad_f, UserDataPtr user_data);

static Bool Eval_G(Index n, Number* x, Bool new_x, Index m, Number* g_, UserDataPtr user_data);

static Bool Eval_Jac_G(Index n, Number *x, Bool new_x,
                              Index m, Index nele_jac,
                              Index *iRow, Index *jCol, Number *values,
                              UserDataPtr user_data);

static Bool Eval_H(Index n, Number *x_, Bool new_x, Number obj_factor,
                          Index m, Number *lambda, Bool new_lambda,
                          Index nele_hess, Index *iRow, Index *jCol,
                          Number *values, UserDataPtr user_data);

static Bool Intermediate(Index alg_mod,
                Index iter_count, Number obj_value,
                Number inf_pr, Number inf_du,
                Number mu, Number d_norm,
                Number regularization_size,
                Number alpha_du, Number alpha_pr,
                Index ls_trials, UserDataPtr user_data);

/******************************/

void getLocalOptimum(datapath dataPath) {

    IpoptProblem problem = setupIpoptProblem();

    Number loglikelihood = INFINITY;

    MyUserData myUserData;
    myUserData.nTheta = getLenTheta();
    myUserData.gradient = alloca(sizeof(double) * myUserData.nTheta);
    myUserData.theta    = alloca(sizeof(double) * myUserData.nTheta);
    myUserData.datapath = dataPath;

    clock_t timeBegin = clock();

    double initialTheta[myUserData.nTheta];
    getInitialTheta(dataPath, initialTheta);

    enum ApplicationReturnStatus status = IpoptSolve(problem, initialTheta, NULL, &loglikelihood, NULL, NULL, NULL, &myUserData);

    clock_t timeEnd = clock();
    double timeElapsed = (double) (timeEnd - timeBegin) / CLOCKS_PER_SEC;

    printf("\n\nIpopt status %d,  final llh: %e, time: %f\n", status, loglikelihood, timeElapsed);

    FreeIpoptProblem(problem);
}

static IpoptProblem setupIpoptProblem()
{
    Index numOptimizationParams = NUM_OPTIMIZATION_PARAMS;

    Number *thetaLowerBounds = alloca(sizeof(Number) * numOptimizationParams);
    fillArray(thetaLowerBounds, numOptimizationParams, 1e-10);

    Number *thetaUpperBounds = alloca(sizeof(Number) * numOptimizationParams);
    fillArray(thetaUpperBounds, numOptimizationParams, 1E4);

    Index numberConstraints = 0;

    Index numNonZeroElementsConstraintJacobian = 0; // TODO only nonzero elements
    Index numNonZeroElementsLagrangianHessian = 0; //NUM_OPTIMIZATION_PARAMS * NUM_OPTIMIZATION_PARAMS; // TODO

    IpoptProblem nlp = CreateIpoptProblem(numOptimizationParams, thetaLowerBounds, thetaUpperBounds,
                                              numberConstraints, NULL, NULL,
                                              numNonZeroElementsConstraintJacobian, numNonZeroElementsLagrangianHessian, 0,
                                              &Eval_F, &Eval_G, &Eval_Grad_F, &Eval_Jac_G, &Eval_H);
    assert(nlp != 0);

    AddIpoptIntOption(nlp, "print_level", 5);
    AddIpoptStrOption(nlp, "print_user_options", "yes");

    AddIpoptStrOption(nlp, "hessian_approximation", "limited-memory");
    AddIpoptStrOption(nlp, "limited_memory_update_type", "bfgs");

    AddIpoptIntOption(nlp, "max_iter", 300);
    AddIpoptNumOption(nlp, "tol", 1e-9);

    AddIpoptIntOption(nlp, "acceptable_iter", 1);
    AddIpoptNumOption(nlp, "acceptable_constr_viol_tol", 1e20);
    AddIpoptNumOption(nlp, "acceptable_dual_inf_tol", 1e20);
    AddIpoptNumOption(nlp, "acceptable_compl_inf_tol", 1e20);
    AddIpoptNumOption(nlp, "acceptable_obj_change_tol", 1e20);

    // TODO check further limited memory options http://www.coin-or.org/Ipopt/documentation/node53.html#opt:hessian_approximation

    OpenIpoptOutputFile(nlp, IPTOPT_LOG_FILE, 6);
    SetIntermediateCallback(nlp, &Intermediate);

    return nlp;
}

/** Type defining the callback function for evaluating the value of
 *  the objective function.  Return value should be set to false if
 *  there was a problem doing the evaluation. */

static Bool Eval_F(Index n, Number *x, Bool new_x, Number *obj_value, UserDataPtr user_data)
{
    static int numFunctionCalls = 0;
    printf("Eval_F (%d) #%d\n", new_x, ++numFunctionCalls);
    fflush(stdout);

    clock_t timeBegin = clock();

    int status = 0;

    MyUserData *data = (MyUserData *) user_data;

    status = evaluateObjectiveFunction(x, n, data->datapath, obj_value, NULL);
    data->objectiveFunctionValue = *obj_value;
    for(int i = 0; i < n; ++i)
        data->theta[i] = x[i];
    fillArray(data->gradient, n, NAN);
    clock_t timeEnd = clock();
    double timeElapsed = (double) (timeEnd - timeBegin) / CLOCKS_PER_SEC;

    logLocalOptimizerObjectiveFunctionEvaluation(data->datapath, numFunctionCalls, x, data->objectiveFunctionValue, timeElapsed, n);

    printf("Done Eval_F (%d) #%d %f\n", new_x, ++numFunctionCalls, *obj_value);

    return !isnan(data->objectiveFunctionValue) && status == 0;
}

/** Type defining the callback function for evaluating the gradient of
 *  the objective function.  Return value should be set to false if
 *  there was a problem doing the evaluation. */

static Bool Eval_Grad_F(Index n, Number *x, Bool new_x, Number *grad_f, UserDataPtr user_data)
{
    static int numFunctionCalls = 0;
    printf("Eval_Grad_F (%d) #%d\n", new_x, ++numFunctionCalls);
    fflush(stdout);

    clock_t timeBegin = clock();

    int status = 0;

    MyUserData *data = (MyUserData *) user_data;

    status = evaluateObjectiveFunction(x, n, data->datapath, &data->objectiveFunctionValue, data->gradient);
    for(int i = 0; i < n; ++i) {
        data->theta[i] = x[i];
        grad_f[i] = data->gradient[i];
    }

    clock_t timeEnd = clock();
    double timeElapsed = (double) (timeEnd - timeBegin) / CLOCKS_PER_SEC;

    logLocalOptimizerObjectiveFunctionGradientEvaluation(data->datapath, numFunctionCalls, x, data->objectiveFunctionValue, data->gradient, timeElapsed, n);


    return !isnan(data->objectiveFunctionValue) && status == 0;
}

/** Type defining the callback function for evaluating the value of
 *  the constraint functions.  Return value should be set to false if
 *  there was a problem doing the evaluation. */

static Bool Eval_G(Index n, Number *x_, Bool new_x, Index m, Number *g_, UserDataPtr user_data)
{
    // no constraints, should never get here
    assert(false);
    return true;
}

/** Type defining the callback function for evaluating the Jacobian of
 *  the constrant functions.  Return value should be set to false if
 *  there was a problem doing the evaluation. */

static Bool Eval_Jac_G(Index n, Number *x, Bool new_x, Index m, Index nele_jac, Index *iRow, Index *jCol, Number *values, UserDataPtr user_data)
{
    // no constraints, nothing to do here, but will be called once

    return true;
}

/** Type defining the callback function for evaluating the Hessian of
 *  the Lagrangian function.  Return value should be set to false if
 *  there was a problem doing the evaluation. */

static Bool Eval_H(Index n, Number *x_, Bool new_x, Number obj_factor, Index m, Number *lambda, Bool new_lambda, Index nele_hess, Index *iRow, Index *jCol, Number *values, UserDataPtr user_data)
{
    static int numFunctionCalls = 0;
    printf("Eval_H #%d\n", ++numFunctionCalls);
    fflush(stdout);

    assert(1==3);
    // TODO not yet used. wait for 2nd order adjoint sensitivities

    assert(m == 0);

    if(iRow && jCol) {
        // provide structure
        // 0 ... nele_hess
    } else {
        // provide values (lower left)
        // 0 ... nele_hess
    }

    // TODO check new_x

    return true;
}

/** Type defining the callback function for giving intermediate
 *  execution control to the user.  If set, it is called once per
 *  iteration, providing the user with some information on the state
 *  of the optimization.  This can be used to print some
 *  user-defined output.  It also gives the user a way to terminate
 *  the optimization prematurely.  If this method returns false,
 *  Ipopt will terminate the optimization. */
/* alg_mod: 0 is regular, 1 is resto */

static Bool Intermediate(Index alg_mod, Index iter_count, Number obj_value, Number inf_pr, Number inf_du, Number mu, Number d_norm, Number regularization_size, Number alpha_du, Number alpha_pr, Index ls_trials, UserDataPtr user_data)
{
    MyUserData *data = (MyUserData *) user_data;

    data->datapath.idxLocalOptimizationIteration = iter_count;

    logLocalOptimizerIteration(data->datapath, iter_count, data->theta, data->objectiveFunctionValue, data->gradient, 0, data->nTheta);

    if(caughtTerminationSignal) {
        fprintf(stderr, "\nCAUGHT SIGTERM... EXITING\n");
        return false;
    }

    return true;
}
