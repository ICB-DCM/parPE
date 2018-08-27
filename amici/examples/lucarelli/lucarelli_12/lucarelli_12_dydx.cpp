#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
using amici::realtype;
#include <cmath> 


#include "species.h"
#include "parameter.h"
#include "fixed_parameter.h"

void dydx_lucarelli_12(double *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h){
    dydx[252] = 1;
    dydx[265] = 1;
    dydx[278] = 1;
    dydx[291] = 1;
    dydx[304] = 1;
    dydx[317] = 1;
    dydx[330] = 1;
    dydx[343] = 1;
    dydx[356] = 1;
    dydx[369] = 1;
    dydx[382] = 1;
    dydx[395] = 1;
}