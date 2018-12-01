#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
using amici::realtype;
#include <cmath>


#include "y.h"
#include "my.h"
#include "p.h"
#include "k.h"
#include "sigmay.h"

void dJydsigmay_model_steadystate_scaled(double *dJydsigmay, const int iy, const realtype *p, const realtype *k, const double *y, const double *sigmay, const double *my){
    switch(iy) {
        case 0:
            dJydsigmay[0] = 1.0/sigmaobservable_x1 - 1.0*pow(-mobservable_x1 + observable_x1, 2)/pow(sigmaobservable_x1, 3);
            break;
        case 1:
            dJydsigmay[1] = 1.0/sigmaobservable_x2 - 1.0*pow(-mobservable_x2 + observable_x2, 2)/pow(sigmaobservable_x2, 3);
            break;
        case 2:
            dJydsigmay[2] = 1.0/sigmaobservable_x3 - 1.0*pow(-mobservable_x3 + observable_x3, 2)/pow(sigmaobservable_x3, 3);
            break;
        case 3:
            dJydsigmay[3] = 1.0/sigmaobservable_x1_scaled - 1.0*pow(-mobservable_x1_scaled + observable_x1_scaled, 2)/pow(sigmaobservable_x1_scaled, 3);
            break;
        case 4:
            dJydsigmay[4] = 1.0/sigmaobservable_x2_offsetted - 1.0*pow(-mobservable_x2_offsetted + observable_x2_offsetted, 2)/pow(sigmaobservable_x2_offsetted, 3);
            break;
        case 5:
            dJydsigmay[5] = 1.0/sigmaobservable_x1withsigma - 1.0*pow(-mobservable_x1withsigma + observable_x1withsigma, 2)/pow(sigmaobservable_x1withsigma, 3);
            break;
}
}