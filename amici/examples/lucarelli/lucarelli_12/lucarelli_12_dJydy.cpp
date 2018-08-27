#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
using amici::realtype;
#include <cmath> 


#include "parameter.h"
#include "fixed_parameter.h"
#include "observable.h"
#include "my.h"
#include "sigmay.h"

void dJydy_lucarelli_12(double *dJydy, const int iy, const realtype *p, const realtype *k, const double *y, const double *sigmay, const double *my){
    switch(iy) {
        case 0:
            dJydy[0] = 1.0*(-mobservable_Ski + observable_Ski)/pow(sigmaobservable_Ski, 2);
            break;
        case 1:
            dJydy[1] = 1.0*(-mobservable_Skil + observable_Skil)/pow(sigmaobservable_Skil, 2);
            break;
        case 2:
            dJydy[2] = 1.0*(-mobservable_Dnmt3a + observable_Dnmt3a)/pow(sigmaobservable_Dnmt3a, 2);
            break;
        case 3:
            dJydy[3] = 1.0*(-mobservable_Sox4 + observable_Sox4)/pow(sigmaobservable_Sox4, 2);
            break;
        case 4:
            dJydy[4] = 1.0*(-mobservable_Jun + observable_Jun)/pow(sigmaobservable_Jun, 2);
            break;
        case 5:
            dJydy[5] = 1.0*(-mobservable_Smad7 + observable_Smad7)/pow(sigmaobservable_Smad7, 2);
            break;
        case 6:
            dJydy[6] = 1.0*(-mobservable_Klf10 + observable_Klf10)/pow(sigmaobservable_Klf10, 2);
            break;
        case 7:
            dJydy[7] = 1.0*(-mobservable_Bmp4 + observable_Bmp4)/pow(sigmaobservable_Bmp4, 2);
            break;
        case 8:
            dJydy[8] = 1.0*(-mobservable_Cxcl15 + observable_Cxcl15)/pow(sigmaobservable_Cxcl15, 2);
            break;
        case 9:
            dJydy[9] = 1.0*(-mobservable_Dusp5 + observable_Dusp5)/pow(sigmaobservable_Dusp5, 2);
            break;
        case 10:
            dJydy[10] = 1.0*(-mobservable_Tgfa + observable_Tgfa)/pow(sigmaobservable_Tgfa, 2);
            break;
        case 11:
            dJydy[11] = 1.0*(-mobservable_Pdk4 + observable_Pdk4)/pow(sigmaobservable_Pdk4, 2);
            break;
}
}