#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
using amici::realtype;
#include <cmath> 


#include "species.h"
#include "parameter.h"
#include "fixed_parameter.h"
#include "flux.h"
#include "dwdx.h"

void J_model_steadystate_scaled(realtype *J, const realtype t, const realtype *x, const double *p, const double *k, const realtype *h, const realtype *w, const realtype *dwdx){
    J[0] = 0.0 - 2.0*dwdx0 - 1.0*dwdx1;
    J[1] = 0.0 - 1.0*dwdx2 + 2.0*dwdx3;
    J[2] = 1.0*dwdx4;
    J[3] = 0.0 + 1.0*dwdx0 - 1.0*dwdx1;
    J[4] = 0.0 - 1.0*dwdx2 - 1.0*dwdx3;
    J[5] = 1.0*dwdx4;
    J[6] = 1.0*dwdx1;
    J[7] = 1.0*dwdx2;
    J[8] = -1.0*dwdx4 - 1.0*dwdx5;
}