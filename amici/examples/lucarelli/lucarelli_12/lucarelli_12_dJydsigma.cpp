#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
using amici::realtype;
#include <cmath> 


#include "parameter.h"
#include "fixed_parameter.h"
#include "observable.h"
#include "my.h"
#include "sigmay.h"

void dJydsigma_lucarelli_12(double *dJydsigma, const int iy, const realtype *p, const realtype *k, const double *y, const double *sigmay, const double *my){
    switch(iy) {
        case 0:
            dJydsigma[0] = -1.0*pow(-mobservable_Ski + observable_Ski, 2)/pow(sigmaobservable_Ski, 3) + 1.0*pow(sigmaobservable_Ski, -1);
            break;
        case 1:
            dJydsigma[1] = -1.0*pow(-mobservable_Skil + observable_Skil, 2)/pow(sigmaobservable_Skil, 3) + 1.0*pow(sigmaobservable_Skil, -1);
            break;
        case 2:
            dJydsigma[2] = -1.0*pow(-mobservable_Dnmt3a + observable_Dnmt3a, 2)/pow(sigmaobservable_Dnmt3a, 3) + 1.0*pow(sigmaobservable_Dnmt3a, -1);
            break;
        case 3:
            dJydsigma[3] = -1.0*pow(-mobservable_Sox4 + observable_Sox4, 2)/pow(sigmaobservable_Sox4, 3) + 1.0*pow(sigmaobservable_Sox4, -1);
            break;
        case 4:
            dJydsigma[4] = -1.0*pow(-mobservable_Jun + observable_Jun, 2)/pow(sigmaobservable_Jun, 3) + 1.0*pow(sigmaobservable_Jun, -1);
            break;
        case 5:
            dJydsigma[5] = -1.0*pow(-mobservable_Smad7 + observable_Smad7, 2)/pow(sigmaobservable_Smad7, 3) + 1.0*pow(sigmaobservable_Smad7, -1);
            break;
        case 6:
            dJydsigma[6] = -1.0*pow(-mobservable_Klf10 + observable_Klf10, 2)/pow(sigmaobservable_Klf10, 3) + 1.0*pow(sigmaobservable_Klf10, -1);
            break;
        case 7:
            dJydsigma[7] = -1.0*pow(-mobservable_Bmp4 + observable_Bmp4, 2)/pow(sigmaobservable_Bmp4, 3) + 1.0*pow(sigmaobservable_Bmp4, -1);
            break;
        case 8:
            dJydsigma[8] = -1.0*pow(-mobservable_Cxcl15 + observable_Cxcl15, 2)/pow(sigmaobservable_Cxcl15, 3) + 1.0*pow(sigmaobservable_Cxcl15, -1);
            break;
        case 9:
            dJydsigma[9] = -1.0*pow(-mobservable_Dusp5 + observable_Dusp5, 2)/pow(sigmaobservable_Dusp5, 3) + 1.0*pow(sigmaobservable_Dusp5, -1);
            break;
        case 10:
            dJydsigma[10] = -1.0*pow(-mobservable_Tgfa + observable_Tgfa, 2)/pow(sigmaobservable_Tgfa, 3) + 1.0*pow(sigmaobservable_Tgfa, -1);
            break;
        case 11:
            dJydsigma[11] = -1.0*pow(-mobservable_Pdk4 + observable_Pdk4, 2)/pow(sigmaobservable_Pdk4, 3) + 1.0*pow(sigmaobservable_Pdk4, -1);
            break;
}
}