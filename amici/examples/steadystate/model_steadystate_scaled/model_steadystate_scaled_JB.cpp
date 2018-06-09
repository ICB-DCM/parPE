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

void JB_model_steadystate_scaled(realtype *JB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx){
    JB[0] = 0.0 - 2.0*dwdx0 - 1.0*dwdx1;
    JB[1] = 0.0 + 1.0*dwdx0 - 1.0*dwdx1;
    JB[2] = 1.0*dwdx1;
    JB[3] = 0.0 - 1.0*dwdx2 + 2.0*dwdx3;
    JB[4] = 0.0 - 1.0*dwdx2 - 1.0*dwdx3;
    JB[5] = 1.0*dwdx2;
    JB[6] = 1.0*dwdx4;
    JB[7] = 1.0*dwdx4;
    JB[8] = -1.0*dwdx4 - 1.0*dwdx5;
}