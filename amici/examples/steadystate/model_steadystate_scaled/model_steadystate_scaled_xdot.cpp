#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
using amici::realtype;
#include <cmath> 


#include "species.h"
#include "parameter.h"
#include "fixed_parameter.h"
#include "flux.h"

void xdot_model_steadystate_scaled(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    xdot[0] = -2.0*w0 - 1.0*w1 + 2.0*w2 + 1.0*w3 + 1.0*w5;
    xdot[1] = 1.0*w0 - 1.0*w1 - 1.0*w2 + 1.0*w3;
    xdot[2] = 1.0*w1 - 1.0*w3 - 1.0*w4;
}