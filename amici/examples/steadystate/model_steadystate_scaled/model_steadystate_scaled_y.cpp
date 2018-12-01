#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
using amici::realtype;
#include <cmath>


#include "w.h"
#include "x.h"
#include "p.h"
#include "k.h"

void y_model_steadystate_scaled(double *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = x1;
    y[1] = x2;
    y[2] = x3;
    y[3] = scaling_x1*x1;
    y[4] = offset_x2 + x2;
    y[5] = x1;
}