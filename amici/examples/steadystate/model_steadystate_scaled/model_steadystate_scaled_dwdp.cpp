#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
using amici::realtype;
#include <cmath> 


#include "species.h"
#include "parameter.h"
#include "fixed_parameter.h"
#include "flux.h"

void dwdp_model_steadystate_scaled(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    dwdp[0] = pow(x1, 2);
    dwdp[1] = x2*x1;
    dwdp[2] = x2;
    dwdp[3] = x3;
    dwdp[4] = 1;
}