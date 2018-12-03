#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
using amici::realtype;
#include <cmath>


#include "w.h"
#include "x.h"
#include "p.h"
#include "k.h"
#include "dwdp.h"
#include "xB.h"

void qBdot_model_steadystate_scaled(realtype *qBdot, const int ip, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdp){
    switch(ip) {
        case 0:
            qBdot[0] = 2.0*dwdp0*xB0 - 1.0*dwdp0*xB1;
            break;
        case 1:
            qBdot[0] = 1.0*dwdp1*xB0 + 1.0*dwdp1*xB1 - 1.0*dwdp1*xB2;
            break;
        case 2:
            qBdot[0] = -2.0*dwdp2*xB0 + 1.0*dwdp2*xB1;
            break;
        case 3:
            qBdot[0] = -1.0*dwdp3*xB0 - 1.0*dwdp3*xB1 + 1.0*dwdp3*xB2;
            break;
        case 4:
            qBdot[0] = -1.0*dwdp4*xB0;
            break;
        case 5:
            break;
        case 6:
            break;
        case 7:
            break;
}
}