#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
using amici::realtype;
#include <cmath>


#include "w.h"
#include "x.h"
#include "p.h"
#include "k.h"
#include "dwdx.h"
#include "JSparse.h"
#include "sx.h"
#include "dxdotdp.h"

void sxdot_model_steadystate_scaled(realtype *sxdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const realtype *sx, const realtype *w, const realtype *dwdx, const realtype *JSparse, const realtype *dxdotdp){
    sxdot[0] = JSparse0*sx0 + JSparse3*sx1 + JSparse6*sx2 + dxdotdp0;
    sxdot[1] = JSparse1*sx0 + JSparse4*sx1 + JSparse7*sx2 + dxdotdp1;
    sxdot[2] = JSparse2*sx0 + JSparse5*sx1 + JSparse8*sx2 + dxdotdp2;
}