#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
using amici::realtype;
#include <cmath> 


#include "species.h"
#include "parameter.h"
#include "fixed_parameter.h"

void w_model_steadystate_scaled(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h){
    w[0] = p1*pow(x1, 2);
    w[1] = p2*x2*x1;
    w[2] = p3*x2;
    w[3] = p4*x3;
    w[4] = k0*x3;
    w[5] = p5;
}