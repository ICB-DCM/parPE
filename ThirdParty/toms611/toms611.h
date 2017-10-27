#ifndef TOMS611_H
#define TOMS611_H

#ifdef __cplusplus
extern "C" {
#endif

#include "f2c.h"
const integer toms611_sumsl_iv_min_length = 60;

constexpr integer toms611_sumsl_v_min_length(integer numOptimizationVariables)
{
    return 71+numOptimizationVariables*(numOptimizationVariables+15)/2;
}

/* see toms611.cpp for details */
enum toms611ReturnCodes {
    x_convergence = 3,
    relative_function_convergence = 4,
    x_and_relative_function_convergence = 5,
    absolute_function_convergence = 6,
    singular_convergence = 7,
    false_convergence = 8,
    max_function_evalulations_exceeded = 9,
    max_iterations_exceeded = 10,
    user_interrupt = 11,
    storage_allocation_done = 14,
    first_error_code = 17,
    attempted_restart_with_changed_variable_number = 17,
    illegal_negative_scaling_factors = 18,
    start_input_value_out_of_range = 19,
    end_input_value_out_of_range = 43,
    infeasible_initial_point = 63,
    bad_parameters = 64,
    gradient_evaluation_failed = 65,
    bad_first_parameter_to_deflt = 67,
    first_input_value_out_of_range = 80,
    invalid_optimization_variable_number = 81
};



/** @brief sumsl_ minimize general unconstrained objective function using analytic gradient and hessian approx. from secant update
 * see toms611.cpp for details
  * @param n in: number of optimization variables
  * @param d in/out: scaling of x
  * @param x  in: initial point; out: optimal point
  * @param calcf cost function: void (*calcf)(int &n, double *x, int &nf, double &f,
       int *uiparm, double *urparm, void *ufparm)
  * @param calcg gradient: void (*calcf)(int &n, double *x, int &nf, double *g,
       int *uiparm, double *urparm, void *ufparm)
  * @param iv in/out: integer options of sumsl
  * @param liv in: length of iv
  * @param lv: length of v
  * @param v: float options of sumsl
  * @param uiparm in/out: user-defined data that will be passed through to calcf/calcg
  * @param urparm
  * @param ufparm
  * @return always 0
  */
int sumsl_(integer &n, doublereal *d, doublereal *x, S_fp
           calcf, S_fp calcg, integer *iv, integer &liv, integer &lv, doublereal
           *v, integer *uiparm, doublereal *urparm, U_fp ufparm);

//int deflt_(integer *alg, integer *iv, integer *liv, integer *lv, doublereal *v);

#ifdef __cplusplus
}
#endif

#endif
