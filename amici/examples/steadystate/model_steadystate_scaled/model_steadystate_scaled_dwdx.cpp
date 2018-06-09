#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
using amici::realtype;
#include <cmath> 


#include "species.h"
#include "parameter.h"
#include "fixed_parameter.h"
#include "flux.h"

void dwdx_model_steadystate_scaled(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    dwdx[0] = 2*p1*x1;
    dwdx[1] = p2*x2;
    dwdx[2] = p2*x1;
    dwdx[3] = p3;
    dwdx[4] = p4;
    dwdx[5] = k0;
}