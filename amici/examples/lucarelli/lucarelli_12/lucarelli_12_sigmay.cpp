#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
using amici::realtype;
#include <cmath> 


#include "parameter.h"
#include "fixed_parameter.h"

void sigmay_lucarelli_12(double *sigmay, const realtype t, const realtype *p, const realtype *k){
    sigmay[0] = 1.0;
    sigmay[1] = 1.0;
    sigmay[2] = 1.0;
    sigmay[3] = 1.0;
    sigmay[4] = 1.0;
    sigmay[5] = 1.0;
    sigmay[6] = 1.0;
    sigmay[7] = 1.0;
    sigmay[8] = 1.0;
    sigmay[9] = 1.0;
    sigmay[10] = 1.0;
    sigmay[11] = 1.0;
}