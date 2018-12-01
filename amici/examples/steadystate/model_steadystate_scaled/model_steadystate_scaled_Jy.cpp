#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
using amici::realtype;
#include <cmath>


#include "y.h"
#include "my.h"
#include "p.h"
#include "k.h"
#include "sigmay.h"

void Jy_model_steadystate_scaled(double *Jy, const int iy, const realtype *p, const realtype *k, const double *y, const double *sigmay, const double *my){
    switch(iy) {
        case 0:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmaobservable_x1, 2)) + 0.5*pow(-mobservable_x1 + observable_x1, 2)/pow(sigmaobservable_x1, 2);
            break;
        case 1:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmaobservable_x2, 2)) + 0.5*pow(-mobservable_x2 + observable_x2, 2)/pow(sigmaobservable_x2, 2);
            break;
        case 2:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmaobservable_x3, 2)) + 0.5*pow(-mobservable_x3 + observable_x3, 2)/pow(sigmaobservable_x3, 2);
            break;
        case 3:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmaobservable_x1_scaled, 2)) + 0.5*pow(-mobservable_x1_scaled + observable_x1_scaled, 2)/pow(sigmaobservable_x1_scaled, 2);
            break;
        case 4:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmaobservable_x2_offsetted, 2)) + 0.5*pow(-mobservable_x2_offsetted + observable_x2_offsetted, 2)/pow(sigmaobservable_x2_offsetted, 2);
            break;
        case 5:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmaobservable_x1withsigma, 2)) + 0.5*pow(-mobservable_x1withsigma + observable_x1withsigma, 2)/pow(sigmaobservable_x1withsigma, 2);
            break;
}
}