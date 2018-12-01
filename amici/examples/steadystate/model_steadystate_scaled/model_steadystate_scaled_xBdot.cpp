#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
using amici::realtype;
#include <cmath>


#include "w.h"
#include "x.h"
#include "p.h"
#include "k.h"
#include "dwdx.h"
#include "xB.h"

void xBdot_model_steadystate_scaled(realtype *xBdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx){
    xBdot[0] = -1.0*dwdx1*xB2 + xB0*(2.0*dwdx0 + 1.0*dwdx1) + xB1*(-1.0*dwdx0 + 1.0*dwdx1);
    xBdot[1] = -1.0*dwdx2*xB2 + xB0*(1.0*dwdx2 - 2.0*dwdx3) + xB1*(1.0*dwdx2 + 1.0*dwdx3);
    xBdot[2] = -1.0*dwdx4*xB0 - 1.0*dwdx4*xB1 + xB2*(1.0*dwdx4 + 1.0*dwdx5);
}