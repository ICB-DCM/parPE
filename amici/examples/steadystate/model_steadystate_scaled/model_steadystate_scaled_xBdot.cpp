#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
using amici::realtype;
#include <cmath> 


#include "species.h"
#include "parameter.h"
#include "fixed_parameter.h"
#include "adjoint.h"
#include "flux.h"
#include "dwdx.h"

void xBdot_model_steadystate_scaled(realtype *xBdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx){
    xBdot[0] = -xB0*(0.0 - 2.0*dwdx0 - 1.0*dwdx1) - xB1*(0.0 + 1.0*dwdx0 - 1.0*dwdx1) - 1.0*xB2*dwdx1;
    xBdot[1] = -xB0*(0.0 - 1.0*dwdx2 + 2.0*dwdx3) - 1.0*xB2*dwdx2 - (0.0 - 1.0*dwdx2 - 1.0*dwdx3)*xB1;
    xBdot[2] = -1.0*xB0*dwdx4 - 1.0*xB1*dwdx4 - (-1.0*dwdx4 - 1.0*dwdx5)*xB2;
}