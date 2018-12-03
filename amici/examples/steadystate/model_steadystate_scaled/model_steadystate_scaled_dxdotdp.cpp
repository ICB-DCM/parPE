#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
using amici::realtype;
#include <cmath>


#include "w.h"
#include "x.h"
#include "p.h"
#include "k.h"
#include "dwdp.h"

void dxdotdp_model_steadystate_scaled(realtype *dxdotdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const realtype *w, const realtype *dwdp){
    switch(ip) {
        case 0:
            dxdotdp[0] = -2.0*dwdp0;
            dxdotdp[1] = 1.0*dwdp0;
            break;
        case 1:
            dxdotdp[0] = -1.0*dwdp1;
            dxdotdp[1] = -1.0*dwdp1;
            dxdotdp[2] = 1.0*dwdp1;
            break;
        case 2:
            dxdotdp[0] = 2.0*dwdp2;
            dxdotdp[1] = -1.0*dwdp2;
            break;
        case 3:
            dxdotdp[0] = 1.0*dwdp3;
            dxdotdp[1] = 1.0*dwdp3;
            dxdotdp[2] = -1.0*dwdp3;
            break;
        case 4:
            dxdotdp[0] = 1.0*dwdp4;
            break;
        case 5:
            break;
        case 6:
            break;
        case 7:
            break;
}
}