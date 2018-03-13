
#include "amici/symbolic_functions.h"
#include "amici/amici_defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

void root_model_neuron_o2(realtype *root, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h) {
  root[0] = x[0]-3.0E1;
}

