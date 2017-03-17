#include "localOptimizationIpopt.h"
#include "optimizationProblem.h"

#include <stdbool.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include <signal.h>
#include <alloca.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include <IpStdCInterface.h>

#include "logging.h"

#define IPTOPT_LOG_FILE "/home/dweindl/src/CanPathProSSH/dw/ipopt.log"

#ifdef INSTALL_SIGNAL_HANDLER
extern volatile sig_atomic_t caughtTerminationSignal;
#endif

static IpoptProblem setupIpoptProblem(OptimizationProblem *problem);

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

int getLocalOptimumIpopt(OptimizationProblem *problem) {

    double *parameters = (double *) malloc(sizeof(*parameters) * problem->numOptimizationParameters);
    // copy, because will be update each iteration
    memcpy(parameters, problem->initialParameters, sizeof(*parameters) * problem->numOptimizationParameters);

    Number loglikelihood = INFINITY;

    IpoptProblem ipoptProblem = setupIpoptProblem(problem);

    clock_t timeBegin = clock();

    enum ApplicationReturnStatus status = IpoptSolve(ipoptProblem,
                                                     parameters,
                                                     NULL,
                                                     &loglikelihood,
                                                     NULL, NULL, NULL,
                                                     problem);

    clock_t timeEnd = clock();
    double timeElapsed = (double) (timeEnd - timeBegin) / CLOCKS_PER_SEC;

    if(problem->logOptimizerFinished)
        problem->logOptimizerFinished(problem, loglikelihood, parameters, timeElapsed, status);

    FreeIpoptProblem(ipoptProblem);
    free(parameters);

    return status < Maximum_Iterations_Exceeded;
}


static IpoptProblem setupIpoptProblem(OptimizationProblem *problem)
{
    assert(sizeof(double) == sizeof(Number));

    // TODO
    Index numberConstraints = 0;
    Index numNonZeroElementsConstraintJacobian = 0; // TODO only nonzero elements
    Index numNonZeroElementsLagrangianHessian = 0; //NUM_OPTIMIZATION_PARAMS * NUM_OPTIMIZATION_PARAMS; // TODO

    IpoptProblem nlp = CreateIpoptProblem(problem->numOptimizationParameters,
                                          problem->parametersMin, problem->parametersMax,
                                              numberConstraints, NULL, NULL,
                                              numNonZeroElementsConstraintJacobian, numNonZeroElementsLagrangianHessian, 0,
                                              &Eval_F, &Eval_G, &Eval_Grad_F, &Eval_Jac_G, &Eval_H);
    assert(nlp != 0);

    if(problem->printToStdout) {
        AddIpoptIntOption(nlp, "print_level", 5);
        AddIpoptStrOption(nlp, "print_user_options", "yes");
    } else {
        AddIpoptIntOption(nlp, "print_level", 0);
        AddIpoptStrOption(nlp, "print_user_options", "no");
        AddIpoptStrOption(nlp, "sb", "yes"); // suppress copyright message
    }

    //    AddIpoptStrOption(nlp, "derivative_test", "first-order");
    //    AddIpoptIntOption(nlp, "derivative_test_first_index", 4130);

    AddIpoptStrOption(nlp, "hessian_approximation", "limited-memory");
    AddIpoptStrOption(nlp, "limited_memory_update_type", "bfgs");

    AddIpoptIntOption(nlp, "max_iter", problem->maxOptimizerIterations);
    AddIpoptNumOption(nlp, "tol", 1e-9);

    //    AddIpoptIntOption(nlp, "acceptable_iter", 1);
    //    AddIpoptNumOption(nlp, "acceptable_constr_viol_tol", 1e20);
    //    AddIpoptNumOption(nlp, "acceptable_dual_inf_tol", 1e20);
    //    AddIpoptNumOption(nlp, "acceptable_compl_inf_tol", 1e20);
    //    AddIpoptNumOption(nlp, "acceptable_obj_change_tol", 1e20);

    // TODO check further limited memory options http://www.coin-or.org/Ipopt/documentation/node53.html#opt:hessian_approximation

    // TODO: remove
    OpenIpoptOutputFile(nlp, IPTOPT_LOG_FILE, 6);
    SetIntermediateCallback(nlp, &Intermediate);

    return nlp;
}

/** Type defining the callback function for evaluating the value of
 *  the objective function.  Return value should be set to false if
 *  there was a problem doing the evaluation. */

static Bool Eval_F(Index n, Number *x, Bool new_x, Number *obj_value, UserDataPtr user_data)
{
    static __thread int numFunctionCalls = 0;
    ++numFunctionCalls;
    // logmessage(LOGLVL_DEBUG, "Eval_F (%d) #%d.", new_x, numFunctionCalls);

    int status = 0;

    clock_t timeBegin = clock();

    OptimizationProblem *problem = (OptimizationProblem *) user_data;
    problem->objectiveFunction(user_data, x, obj_value);

    clock_t timeEnd = clock();
    double timeElapsed = (double) (timeEnd - timeBegin) / CLOCKS_PER_SEC;

    if(problem->logObjectiveFunctionEvaluation)
        problem->logObjectiveFunctionEvaluation(problem, x, *obj_value, numFunctionCalls, timeElapsed);

    return !isnan(*obj_value) && status == 0;
}

/** Type defining the callback function for evaluating the gradient of
 *  the objective function.  Return value should be set to false if
 *  there was a problem doing the evaluation. */

static Bool Eval_Grad_F(Index n, Number *x, Bool new_x, Number *grad_f, UserDataPtr user_data)
{
    static __thread int numFunctionCalls = 0;
    ++numFunctionCalls;
    // logmessage(LOGLVL_DEBUG, "Eval_Grad_F (%d) #%d", new_x, numFunctionCalls);

    int status = 0;

    clock_t timeBegin = clock();

    OptimizationProblem *problem = (OptimizationProblem *) user_data;
    double objectiveFunctionValue;
    problem->objectiveFunctionGradient(user_data, x, &objectiveFunctionValue, grad_f);

    clock_t timeEnd = clock();
    double timeElapsed = (double) (timeEnd - timeBegin) / CLOCKS_PER_SEC;

    if(problem->logObjectiveFunctionGradientEvaluation)
        problem->logObjectiveFunctionGradientEvaluation(problem, x, objectiveFunctionValue, grad_f, numFunctionCalls, timeElapsed);

    return !isnan(objectiveFunctionValue) && status == 0;
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
    static __thread int numFunctionCalls = 0;
    logmessage(LOGLVL_DEBUG, "Eval_H #%d", ++numFunctionCalls);

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

static Bool Intermediate(Index alg_mod,
                         Index iter_count,
                         Number obj_value,
                         Number inf_pr, Number inf_du,
                         Number mu,
                         Number d_norm,
                         Number regularization_size,
                         Number alpha_du, Number alpha_pr,
                         Index ls_trials,
                         UserDataPtr user_data)
{
    OptimizationProblem *problem = (OptimizationProblem *) user_data;

    int status = true;

    if(problem->intermediateFunction)
            status = problem->intermediateFunction(problem, alg_mod,
                                               iter_count, obj_value,
                                               inf_pr,  inf_du,
                                               mu, d_norm,
                                               regularization_size,
                                               alpha_du,  alpha_pr,
                                               ls_trials);

#ifdef INSTALL_SIGNAL_HANDLER
    if(caughtTerminationSignal) {
        logmessage(LOGLVL_CRITICAL, "CAUGHT SIGTERM... EXITING.");
        return false;
    }
#endif

    return status;
}
