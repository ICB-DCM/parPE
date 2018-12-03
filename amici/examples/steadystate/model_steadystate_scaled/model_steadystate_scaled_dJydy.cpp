#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
using amici::realtype;
#include <cmath>


#include "y.h"
#include "my.h"
#include "p.h"
#include "k.h"
#include "sigmay.h"

void dJydy_model_steadystate_scaled(double *dJydy, const int iy, const realtype *p, const realtype *k, const double *y, const double *sigmay, const double *my){
    switch(iy) {
        case 0:
            dJydy[0] = 0.5*(-2*mobservable_x1 + 2*observable_x1)/pow(sigmaobservable_x1, 2);
            break;
        case 1:
            dJydy[1] = 0.5*(-2*mobservable_x2 + 2*observable_x2)/pow(sigmaobservable_x2, 2);
            break;
        case 2:
            dJydy[2] = 0.5*(-2*mobservable_x3 + 2*observable_x3)/pow(sigmaobservable_x3, 2);
            break;
        case 3:
            dJydy[3] = 0.5*(-2*mobservable_x1_scaled + 2*observable_x1_scaled)/pow(sigmaobservable_x1_scaled, 2);
            break;
        case 4:
            dJydy[4] = 0.5*(-2*mobservable_x2_offsetted + 2*observable_x2_offsetted)/pow(sigmaobservable_x2_offsetted, 2);
            break;
        case 5:
            dJydy[5] = 0.5*(-2*mobservable_x1withsigma + 2*observable_x1withsigma)/pow(sigmaobservable_x1withsigma, 2);
            break;
}
}