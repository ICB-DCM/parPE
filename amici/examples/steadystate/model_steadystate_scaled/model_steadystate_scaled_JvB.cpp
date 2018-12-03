#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
using amici::realtype;
#include <cmath>


#include "w.h"
#include "x.h"
#include "p.h"
#include "k.h"
#include "dwdx.h"
#include "vB.h"

void JvB_model_steadystate_scaled(realtype *JvB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *vB, const realtype *w, const realtype *dwdx){
    JvB[0] = 1.0*dwdx1*vB2 + vB0*(-2.0*dwdx0 - 1.0*dwdx1) + vB1*(1.0*dwdx0 - 1.0*dwdx1);
    JvB[1] = 1.0*dwdx2*vB2 + vB0*(-1.0*dwdx2 + 2.0*dwdx3) + vB1*(-1.0*dwdx2 - 1.0*dwdx3);
    JvB[2] = 1.0*dwdx4*vB0 + 1.0*dwdx4*vB1 + vB2*(-1.0*dwdx4 - 1.0*dwdx5);
}