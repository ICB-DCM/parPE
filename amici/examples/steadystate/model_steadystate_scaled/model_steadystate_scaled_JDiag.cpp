#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
using amici::realtype;
#include <cmath> 


#include "species.h"
#include "parameter.h"
#include "fixed_parameter.h"
#include "flux.h"
#include "dwdx.h"

void JDiag_model_steadystate_scaled(realtype *JDiag, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    JDiag[0] = 0.0 - 2.0*dwdx0 - 1.0*dwdx1;
    JDiag[1] = 0.0 - 1.0*dwdx2 - 1.0*dwdx3;
    JDiag[2] = -1.0*dwdx4 - 1.0*dwdx5;
}