#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
using amici::realtype;
#include <cmath>


#include "w.h"
#include "x.h"
#include "p.h"
#include "k.h"

void xdot_model_steadystate_scaled(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    xdot[0] = -2.0*flux_r0 - 1.0*flux_r1 + 2.0*flux_r2 + 1.0*flux_r3 + 1.0*flux_r5;
    xdot[1] = 1.0*flux_r0 - 1.0*flux_r1 - 1.0*flux_r2 + 1.0*flux_r3;
    xdot[2] = 1.0*flux_r1 - 1.0*flux_r3 - 1.0*flux_r4;
}