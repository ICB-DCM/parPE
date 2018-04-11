
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

void dydx_model_steadystate_scaled(double *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h) {
  dydx[0+0*5] = 1.0;
  dydx[1+1*5] = 1.0;
  dydx[2+2*5] = 1.0;
  dydx[3+0*5] = p[6];
  dydx[4+1*5] = 1.0;
}

