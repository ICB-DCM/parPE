#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
using amici::realtype;
#include <cmath> 


#include "species.h"
#include "parameter.h"
#include "fixed_parameter.h"
#include "adjoint.h"
#include "flux.h"
#include "dwdp.h"

void qBdot_model_steadystate_scaled(realtype *qBdot, const int ip, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdp){
    switch(ip) {
        case 0:
            qBdot[0] = 2.0*xB0*dwdp0 - 1.0*xB1*dwdp0;
            break;
        case 1:
            qBdot[0] = 1.0*xB0*dwdp1 + 1.0*xB1*dwdp1 - 1.0*xB2*dwdp1;
            break;
        case 2:
            qBdot[0] = -2.0*xB0*dwdp2 + 1.0*xB1*dwdp2;
            break;
        case 3:
            qBdot[0] = -1.0*xB0*dwdp3 - 1.0*xB1*dwdp3 + 1.0*xB2*dwdp3;
            break;
        case 4:
            qBdot[0] = -1.0*xB0*dwdp4;
            break;
        case 5:
            break;
        case 6:
            break;
}
}