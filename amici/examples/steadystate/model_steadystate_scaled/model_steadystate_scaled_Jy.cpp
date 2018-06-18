#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
using amici::realtype;
#include <cmath> 


#include "parameter.h"
#include "fixed_parameter.h"
#include "observable.h"
#include "my.h"
#include "sigmay.h"

void Jy_model_steadystate_scaled(double *nllh, const int iy, const realtype *p, const realtype *k, const double *y, const double *sigmay, const double *my){
    switch(iy) {
        case 0:
            nllh[0] = 0.5*pow(-mobservable_x1 + observable_x1, 2)/pow(sigmaobservable_x1, 2) + 0.5*log(2*M_PI*pow(sigmaobservable_x1, 2));
            break;
        case 1:
            nllh[0] = 0.5*pow(-mobservable_x2 + observable_x2, 2)/pow(sigmaobservable_x2, 2) + 0.5*log(2*M_PI*pow(sigmaobservable_x2, 2));
            break;
        case 2:
            nllh[0] = 0.5*pow(-mobservable_x3 + observable_x3, 2)/pow(sigmaobservable_x3, 2) + 0.5*log(2*M_PI*pow(sigmaobservable_x3, 2));
            break;
        case 3:
            nllh[0] = 0.5*pow(-mobservable_x1_scaled + observable_x1_scaled, 2)/pow(sigmaobservable_x1_scaled, 2) + 0.5*log(2*M_PI*pow(sigmaobservable_x1_scaled, 2));
            break;
        case 4:
            nllh[0] = 0.5*pow(-mobservable_x2_offsetted + observable_x2_offsetted, 2)/pow(sigmaobservable_x2_offsetted, 2) + 0.5*log(2*M_PI*pow(sigmaobservable_x2_offsetted, 2));
            break;
        case 5:
            nllh[0] = 0.5*pow(-mobservable_x1withsigma + observable_x1withsigma, 2)/pow(sigmaobservable_x1withsigma, 2) + 0.5*log(2*M_PI*pow(sigmaobservable_x1withsigma, 2));
            break;
}
}