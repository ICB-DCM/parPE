#include "localoptimization.h"

#include <stdbool.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include "objectivefunction.h"

void getLocalOptimum(double *initialTheta, loggerdata *datalogger) {

    IpoptProblem problem = setupIpoptProblem();

    Number loglikelihood = INFINITY;

    MyUserData myUserData;
    myUserData.nTheta = NUM_OPTIMIZATION_PARAMS;
    myUserData.gradient = malloc(sizeof(double) * NUM_OPTIMIZATION_PARAMS);
    myUserData.theta    = malloc(sizeof(double) * NUM_OPTIMIZATION_PARAMS);
    myUserData.datalogger = datalogger;

    clock_t timeBegin = clock();

    enum ApplicationReturnStatus status = IpoptSolve(problem, initialTheta, NULL, &loglikelihood, NULL, NULL, NULL, &myUserData);

    free(myUserData.gradient);
    free(myUserData.theta);

    clock_t timeEnd = clock();
    double timeElapsed = (double) (timeEnd - timeBegin) / CLOCKS_PER_SEC;

    printf("Ipopt status %d,  final llh: %e, time: %f\n", status, loglikelihood, timeElapsed);

    FreeIpoptProblem(problem);

    // TODO return solution in initialTheta
    printf("Theta 1..10: ");
    printfArray(initialTheta, 10, "%e ");
    printf("\n");
}

IpoptProblem setupIpoptProblem()
{
    Index numOptimizationParams = NUM_OPTIMIZATION_PARAMS;

    Number *thetaLowerBounds = malloc(sizeof(Number) * numOptimizationParams);
    fillArray(thetaLowerBounds, numOptimizationParams, 1e-10);

    Number *thetaUpperBounds = malloc(sizeof(Number) * numOptimizationParams);
    fillArray(thetaUpperBounds, numOptimizationParams, 1E4);

    Index numberConstraints = 0;

    Index numNonZeroElementsConstraintJacobian = 0; // TODO only nonzero elements
    Index numNonZeroElementsLagrangianHessian = 0; //NUM_OPTIMIZATION_PARAMS * NUM_OPTIMIZATION_PARAMS; // TODO

    IpoptProblem nlp = CreateIpoptProblem(numOptimizationParams, thetaLowerBounds, thetaUpperBounds,
                                              numberConstraints, NULL, NULL,
                                              numNonZeroElementsConstraintJacobian, numNonZeroElementsLagrangianHessian, 0,
                                              &Eval_F, &Eval_G, &Eval_Grad_F, &Eval_Jac_G, &Eval_H);
    assert(nlp != 0);

    free(thetaLowerBounds);
    free(thetaUpperBounds);

    AddIpoptIntOption(nlp, "print_level", 5);
    AddIpoptStrOption(nlp, "print_user_options", "yes");
    AddIpoptNumOption(nlp, "tol", 1e-9);
    AddIpoptIntOption(nlp, "max_iter", 1);
    AddIpoptStrOption(nlp, "hessian_approximation", "limited-memory");
    AddIpoptStrOption(nlp, "limited_memory_update_type", "bfgs");

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

Bool Eval_F(Index n, Number *x, Bool new_x, Number *obj_value, UserDataPtr user_data)
{
    assert(n == NUM_OPTIMIZATION_PARAMS);

    static int numFunctionCalls = 0;
    printf("Eval_F (%d) #%d\n", new_x, ++numFunctionCalls);
    fflush(stdout);

    clock_t timeBegin = clock();

    int status = 0;

    MyUserData *data = (MyUserData *) user_data;

    if(new_x) {
        status = evaluateObjectiveFunction(x, n, NUM_CELL_LINES, obj_value, data->gradient);
        data->objectiveFunctionValue = obj_value[0];
        for(int i = 0; i < n; ++i)
            data->theta[i] = x[i];
    } else {
        obj_value[0] = data->objectiveFunctionValue;
    }

    clock_t timeEnd = clock();
    double timeElapsed = (double) (timeEnd - timeBegin) / CLOCKS_PER_SEC;

    if(data->datalogger) {
        logLocalOptimizerObjectiveFunctionEvaluation(*data->datalogger, numFunctionCalls, x, data->objectiveFunctionValue, data->gradient, timeElapsed, n);
    }

    return status == 0;
}

Bool Eval_Grad_F(Index n, Number *x, Bool new_x, Number *grad_f, UserDataPtr user_data)
{
    assert(n == NUM_OPTIMIZATION_PARAMS);

    static int numFunctionCalls = 0;
    printf("Eval_Grad_F (%d) #%d\n", new_x, ++numFunctionCalls);
    fflush(stdout);

    if(new_x) {
        Number objVal = 0;
        Bool fstatus = Eval_F(n, x, new_x, &objVal, user_data);
        if(!fstatus)
            return FALSE;
    } else {
        MyUserData *data = (MyUserData*) user_data;
        for(int i = 0; i < n; ++i) {
            grad_f[i] = data->gradient[i];
        }
    }
    printfArray(grad_f, n, "%f ");

    return true;
}

Bool Eval_G(Index n, Number *x_, Bool new_x, Index m, Number *g_, UserDataPtr user_data)
{
    // no constraints, should never get here
    assert(false);
}

Bool Eval_Jac_G(Index n, Number *x, Bool new_x, Index m, Index nele_jac, Index *iRow, Index *jCol, Number *values, UserDataPtr user_data)
{
    // no constraints, nothing to do here, but will be called once

    return true;
}

Bool Eval_H(Index n, Number *x_, Bool new_x, Number obj_factor, Index m, Number *lambda, Bool new_lambda, Index nele_hess, Index *iRow, Index *jCol, Number *values, UserDataPtr user_data)
{
    static int numFunctionCalls = 0;
    printf("Eval_H #%d\n", ++numFunctionCalls);
    fflush(stdout);

    assert(1==3);
    // TODO not yet used. wait for 2nd order adjoint sensitivities

    assert(n == NUM_OPTIMIZATION_PARAMS);
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

Bool Intermediate(Index alg_mod, Index iter_count, Number obj_value, Number inf_pr, Number inf_du, Number mu, Number d_norm, Number regularization_size, Number alpha_du, Number alpha_pr, Index ls_trials, UserDataPtr user_data)
{
    MyUserData *data = (MyUserData *) user_data;

    if(data->datalogger) {
        logLocalOptimizerIteration(*data->datalogger, iter_count, data->theta, data->objectiveFunctionValue, data->gradient, 0, data->nTheta);
    }

    return true;
}

