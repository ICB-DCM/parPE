#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
using amici::realtype;
#include <cmath> 


#include "species.h"
#include "sensitivity.h"
#include "parameter.h"
#include "fixed_parameter.h"
#include "flux.h"
#include "dxdotdp.h"
#include "dwdx.h"
#include "JSparse.h"

void sxdot_model_steadystate_scaled(realtype *sxdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const realtype *sx, const realtype *w, const realtype *dwdx, const realtype *J, const realtype *dxdotdp){
    sxdot[0] = dxdotdp0 + J0*sx0 + J3*sx1 + J6*sx2;
    sxdot[1] = dxdotdp1 + J1*sx0 + J4*sx1 + J7*sx2;
    sxdot[2] = dxdotdp2 + J2*sx0 + J5*sx1 + J8*sx2;
}