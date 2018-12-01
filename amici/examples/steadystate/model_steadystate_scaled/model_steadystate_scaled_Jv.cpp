#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
using amici::realtype;
#include <cmath>


#include "w.h"
#include "x.h"
#include "p.h"
#include "k.h"
#include "dwdx.h"
#include "v.h"

void Jv_model_steadystate_scaled(realtype *Jv, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *v, const realtype *w, const realtype *dwdx){
    Jv[0] = 1.0*dwdx4*v2 + v0*(-2.0*dwdx0 - 1.0*dwdx1) + v1*(-1.0*dwdx2 + 2.0*dwdx3);
    Jv[1] = 1.0*dwdx4*v2 + v0*(1.0*dwdx0 - 1.0*dwdx1) + v1*(-1.0*dwdx2 - 1.0*dwdx3);
    Jv[2] = 1.0*dwdx1*v0 + 1.0*dwdx2*v1 + v2*(-1.0*dwdx4 - 1.0*dwdx5);
}