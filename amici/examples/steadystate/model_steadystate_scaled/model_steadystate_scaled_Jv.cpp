#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
using amici::realtype;
#include <cmath> 


#include "species.h"
#include "vector.h"
#include "parameter.h"
#include "fixed_parameter.h"
#include "flux.h"
#include "dwdx.h"

void Jv_model_steadystate_scaled(realtype *Jv, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *v, const realtype *w, const realtype *dwdx){
    Jv[0] = v0*(0.0 - 2.0*dwdx0 - 1.0*dwdx1) + v1*(0.0 - 1.0*dwdx2 + 2.0*dwdx3) + 1.0*v2*dwdx4;
    Jv[1] = v0*(0.0 + 1.0*dwdx0 - 1.0*dwdx1) + 1.0*v2*dwdx4 + (0.0 - 1.0*dwdx2 - 1.0*dwdx3)*v1;
    Jv[2] = 1.0*v0*dwdx1 + 1.0*v1*dwdx2 + (-1.0*dwdx4 - 1.0*dwdx5)*v2;
}