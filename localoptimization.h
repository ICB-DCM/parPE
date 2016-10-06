#ifndef LOCAL_OPTIMIZATION_H
#define LOCAL_OPTIMIZATION_H

#include <IpStdCInterface.h>

#include "dataprovider.h"
#include <include/udata.h>
#include <include/edata.h>
#include <include/rdata.h>

#include "logger.h"

#define IPTOPT_LOG_FILE "/home/dweindl/src/CanPathProSSH/dw/ipopt.log"

void getLocalOptimum(double *initialTheta, loggerdata *datalogger);
IpoptProblem setupIpoptProblem();
void getLocalOptimumOld();

typedef struct {
    int nTheta;
    double *theta;
    UserData udata;
    ExpData edata;
    ReturnData rdata;
    double *gradient;
    double objectiveFunctionValue;
    loggerdata *datalogger;
} MyUserData;

/** Type defining the callback function for evaluating the value of
 *  the objective function.  Return value should be set to false if
 *  there was a problem doing the evaluation. */

Bool Eval_F(Index n, Number* x, Bool new_x, Number* obj_value, UserDataPtr user_data);

/** Type defining the callback function for evaluating the gradient of
 *  the objective function.  Return value should be set to false if
 *  there was a problem doing the evaluation. */

Bool Eval_Grad_F(Index n, Number* x, Bool new_x, Number* grad_f, UserDataPtr user_data);

/** Type defining the callback function for evaluating the value of
 *  the constraint functions.  Return value should be set to false if
 *  there was a problem doing the evaluation. */

Bool Eval_G(Index n, Number* x, Bool new_x, Index m, Number* g_, UserDataPtr user_data);

/** Type defining the callback function for evaluating the Jacobian of
 *  the constrant functions.  Return value should be set to false if
 *  there was a problem doing the evaluation. */

Bool Eval_Jac_G(Index n, Number *x, Bool new_x,
                              Index m, Index nele_jac,
                              Index *iRow, Index *jCol, Number *values,
                              UserDataPtr user_data);

/** Type defining the callback function for evaluating the Hessian of
 *  the Lagrangian function.  Return value should be set to false if
 *  there was a problem doing the evaluation. */
Bool Eval_H(Index n, Number *x_, Bool new_x, Number obj_factor,
                          Index m, Number *lambda, Bool new_lambda,
                          Index nele_hess, Index *iRow, Index *jCol,
                          Number *values, UserDataPtr user_data);

/** Type defining the callback function for giving intermediate
 *  execution control to the user.  If set, it is called once per
 *  iteration, providing the user with some information on the state
 *  of the optimization.  This can be used to print some
 *  user-defined output.  It also gives the user a way to terminate
 *  the optimization prematurely.  If this method returns false,
 *  Ipopt will terminate the optimization. */
Bool Intermediate(Index alg_mod, /* 0 is regular, 1 is resto */
                Index iter_count, Number obj_value,
                Number inf_pr, Number inf_du,
                Number mu, Number d_norm,
                Number regularization_size,
                Number alpha_du, Number alpha_pr,
                Index ls_trials, UserDataPtr user_data);

#endif
