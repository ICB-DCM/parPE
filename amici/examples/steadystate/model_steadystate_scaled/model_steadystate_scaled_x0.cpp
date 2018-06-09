#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
using amici::realtype;
#include <cmath> 


#include "parameter.h"
#include "fixed_parameter.h"

void x0_model_steadystate_scaled(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[0] = 0.1;
    x0[1] = 0.4;
    x0[2] = 0.7;
}