
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

void x0_model_steadystate_scaled(realtype *x0, const realtype t, const realtype *p, const realtype *k) {
  x0[0] = 1.0/1.0E1;
  x0[1] = 2.0/5.0;
  x0[2] = 7.0/1.0E1;
}

