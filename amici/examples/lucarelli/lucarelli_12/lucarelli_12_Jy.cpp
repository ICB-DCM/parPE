#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
using amici::realtype;
#include <cmath> 


#include "parameter.h"
#include "fixed_parameter.h"
#include "observable.h"
#include "my.h"
#include "sigmay.h"

void Jy_lucarelli_12(double *nllh, const int iy, const realtype *p, const realtype *k, const double *y, const double *sigmay, const double *my){
    switch(iy) {
        case 0:
            nllh[0] = 0.5*pow(-mobservable_Ski + observable_Ski, 2)/pow(sigmaobservable_Ski, 2) + 0.5*log(2*M_PI*pow(sigmaobservable_Ski, 2));
            break;
        case 1:
            nllh[0] = 0.5*pow(-mobservable_Skil + observable_Skil, 2)/pow(sigmaobservable_Skil, 2) + 0.5*log(2*M_PI*pow(sigmaobservable_Skil, 2));
            break;
        case 2:
            nllh[0] = 0.5*pow(-mobservable_Dnmt3a + observable_Dnmt3a, 2)/pow(sigmaobservable_Dnmt3a, 2) + 0.5*log(2*M_PI*pow(sigmaobservable_Dnmt3a, 2));
            break;
        case 3:
            nllh[0] = 0.5*pow(-mobservable_Sox4 + observable_Sox4, 2)/pow(sigmaobservable_Sox4, 2) + 0.5*log(2*M_PI*pow(sigmaobservable_Sox4, 2));
            break;
        case 4:
            nllh[0] = 0.5*pow(-mobservable_Jun + observable_Jun, 2)/pow(sigmaobservable_Jun, 2) + 0.5*log(2*M_PI*pow(sigmaobservable_Jun, 2));
            break;
        case 5:
            nllh[0] = 0.5*pow(-mobservable_Smad7 + observable_Smad7, 2)/pow(sigmaobservable_Smad7, 2) + 0.5*log(2*M_PI*pow(sigmaobservable_Smad7, 2));
            break;
        case 6:
            nllh[0] = 0.5*pow(-mobservable_Klf10 + observable_Klf10, 2)/pow(sigmaobservable_Klf10, 2) + 0.5*log(2*M_PI*pow(sigmaobservable_Klf10, 2));
            break;
        case 7:
            nllh[0] = 0.5*pow(-mobservable_Bmp4 + observable_Bmp4, 2)/pow(sigmaobservable_Bmp4, 2) + 0.5*log(2*M_PI*pow(sigmaobservable_Bmp4, 2));
            break;
        case 8:
            nllh[0] = 0.5*pow(-mobservable_Cxcl15 + observable_Cxcl15, 2)/pow(sigmaobservable_Cxcl15, 2) + 0.5*log(2*M_PI*pow(sigmaobservable_Cxcl15, 2));
            break;
        case 9:
            nllh[0] = 0.5*pow(-mobservable_Dusp5 + observable_Dusp5, 2)/pow(sigmaobservable_Dusp5, 2) + 0.5*log(2*M_PI*pow(sigmaobservable_Dusp5, 2));
            break;
        case 10:
            nllh[0] = 0.5*pow(-mobservable_Tgfa + observable_Tgfa, 2)/pow(sigmaobservable_Tgfa, 2) + 0.5*log(2*M_PI*pow(sigmaobservable_Tgfa, 2));
            break;
        case 11:
            nllh[0] = 0.5*pow(-mobservable_Pdk4 + observable_Pdk4, 2)/pow(sigmaobservable_Pdk4, 2) + 0.5*log(2*M_PI*pow(sigmaobservable_Pdk4, 2));
            break;
}
}