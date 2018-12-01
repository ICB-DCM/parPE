#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
using amici::realtype;
#include <cmath>


#include "w.h"
#include "x.h"
#include "p.h"
#include "k.h"
#include "dwdx.h"

void dydx_model_steadystate_scaled(double *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    dydx[0] = 1;
    dydx[3] = scaling_x1;
    dydx[5] = 1;
    dydx[7] = 1;
    dydx[10] = 1;
    dydx[14] = 1;
}