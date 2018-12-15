#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
using amici::realtype;
#include <cmath>


#include "p.h"
#include "k.h"

void x0_model_steadystate_scaled(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[0] = 0.10000000000000001;
    x0[1] = 0.40000000000000002;
    x0[2] = 0.69999999999999996;
}