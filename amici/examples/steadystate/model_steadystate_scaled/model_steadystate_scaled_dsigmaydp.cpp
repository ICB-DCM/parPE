#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
using amici::realtype;
#include <cmath> 


#include "parameter.h"
#include "fixed_parameter.h"

void dsigmaydp_model_steadystate_scaled(double *dsigmaydp, const realtype t, const realtype *p, const realtype *k, const int ip){
    switch(ip) {
        case 0:
            break;
        case 1:
            break;
        case 2:
            break;
        case 3:
            break;
        case 4:
            break;
        case 5:
            break;
        case 6:
            break;
        case 7:
            dsigmaydp[5] = 1;
            break;
}
}