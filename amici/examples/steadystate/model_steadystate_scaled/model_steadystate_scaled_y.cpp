
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

void y_model_steadystate_scaled(double *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h) {
  y[0] = x[0];
  y[1] = x[1];
  y[2] = x[2];
  y[3] = p[6]*x[0];
  y[4] = p[7]+x[1];
}

