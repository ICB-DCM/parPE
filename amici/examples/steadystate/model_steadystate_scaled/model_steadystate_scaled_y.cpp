#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
using amici::realtype;
#include <cmath> 


#include "species.h"
#include "parameter.h"
#include "fixed_parameter.h"

void y_model_steadystate_scaled(double *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h){
    y[0] = x1;
    y[1] = x2;
    y[2] = x3;
    y[3] = x1*scaling_x1;
    y[4] = offset_x2 + x2;
}