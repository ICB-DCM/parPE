#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
using amici::realtype;
#include <cmath> 


#include "species.h"
#include "vectorB.h"
#include "parameter.h"
#include "fixed_parameter.h"
#include "adjoint.h"
#include "flux.h"
#include "dwdx.h"

void JvB_model_steadystate_scaled(realtype *JvB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *vB, const realtype *w, const realtype *dwdx){
    JvB[0] = vB0*(0.0 - 2.0*dwdx0 - 1.0*dwdx1) + vB1*(0.0 + 1.0*dwdx0 - 1.0*dwdx1) + 1.0*vB2*dwdx1;
    JvB[1] = vB0*(0.0 - 1.0*dwdx2 + 2.0*dwdx3) + 1.0*vB2*dwdx2 + (0.0 - 1.0*dwdx2 - 1.0*dwdx3)*vB1;
    JvB[2] = 1.0*vB0*dwdx4 + 1.0*vB1*dwdx4 + (-1.0*dwdx4 - 1.0*dwdx5)*vB2;
}