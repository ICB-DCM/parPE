#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
using amici::realtype;
#include <cmath> 


#include "species.h"
#include "parameter.h"
#include "fixed_parameter.h"

void y_lucarelli_12(double *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h){
    y[0] = geneA;
    y[1] = geneB;
    y[2] = geneC;
    y[3] = geneD;
    y[4] = geneE;
    y[5] = geneF;
    y[6] = geneG;
    y[7] = geneH;
    y[8] = geneI;
    y[9] = geneJ;
    y[10] = geneK;
    y[11] = geneL;
}