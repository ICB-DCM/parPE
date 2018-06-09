#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
using amici::realtype;
#include <cmath> 


#include "parameter.h"
#include "fixed_parameter.h"

void sigma_y_model_steadystate_scaled(double *sigmay, const realtype t, const realtype *p, const realtype *k){
    sigmay[0] = 1;
    sigmay[1] = 1;
    sigmay[2] = 1;
    sigmay[3] = 1;
    sigmay[4] = 1;
}