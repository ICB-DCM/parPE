
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

void dydp_model_steadystate_scaled(double *dydp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip) {
switch (ip) {
  case 6: {
  dydp[3] = x[0];

  } break;

  case 7: {
  dydp[4] = 1.0;

  } break;

}
}

