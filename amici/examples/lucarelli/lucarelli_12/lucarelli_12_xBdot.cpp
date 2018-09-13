#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
using amici::realtype;
#include <cmath> 


#include "species.h"
#include "parameter.h"
#include "fixed_parameter.h"
#include "adjoint.h"
#include "flux.h"
#include "dwdx.h"

void xBdot_lucarelli_12(realtype *xBdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx){
    xBdot[0] = 0.0759301442672741*xB0*dwdx0 + 0.0759301442672741*xB1*dwdx0 - 0.0759301442672741*xB2*dwdx0;
    xBdot[1] = 0.0759301442672741*xB0*dwdx1 + 0.0759301442672741*xB1*dwdx1 - 0.0759301442672741*xB2*dwdx1;
    xBdot[2] = -0.0759301442672741*xB12*dwdx3 - 0.0759301442672741*xB13*dwdx4 + 0.0759301442672741*xB2*dwdx2 + 0.0759301442672741*xB3*dwdx3 + 0.0759301442672741*xB4*dwdx4;
    xBdot[3] = -0.0759301442672741*xB12*dwdx6 + 0.151860288534548*xB5*dwdx5 - 0.0759301442672741*xB6*dwdx5 - (0.0 - 0.0759301442672741*dwdx5 - 0.0759301442672741*dwdx6)*xB3;
    xBdot[4] = -0.0759301442672741*xB13*dwdx7 + 0.0759301442672741*xB4*dwdx7;
    xBdot[5] = -xB12*(0.0 - 0.151860288534548*dwdx10 - 0.0759301442672741*dwdx12 - 0.0759301442672741*dwdx14) - xB13*(0.0 - 0.151860288534548*dwdx11 - 0.0759301442672741*dwdx13 - 0.0759301442672741*dwdx14) - 0.0759301442672741*xB14*dwdx10 - 0.0759301442672741*xB17*dwdx11 - 0.0759301442672741*xB18*dwdx14 - 0.0759301442672741*xB19*dwdx13 - 0.0759301442672741*xB20*dwdx12 + 0.0759301442672741*xB3*dwdx8 - xB5*(0.0 - 0.0759301442672741*dwdx10 - 0.0759301442672741*dwdx11 - 0.151860288534548*dwdx12 - 0.151860288534548*dwdx13 - 0.0759301442672741*dwdx14 - 0.151860288534548*dwdx8 - 0.227790432801822*dwdx9) - 0.0759301442672741*xB6*dwdx8 - 0.0759301442672741*xB9*dwdx9;
    xBdot[6] = -0.0759301442672741*xB3*dwdx15 - 0.151860288534548*xB5*dwdx15 + 0.0759301442672741*xB6*dwdx15;
    xBdot[7] = -0.0759301442672741*xB10*dwdx16 - 0.151860288534548*xB12*dwdx16 + 0.0759301442672741*xB7*dwdx16;
    xBdot[8] = -0.0759301442672741*xB11*dwdx17 - 0.151860288534548*xB13*dwdx17 + 0.0759301442672741*xB8*dwdx17;
    xBdot[9] = -0.227790432801822*xB5*dwdx18 + 0.0759301442672741*xB9*dwdx18;
    xBdot[10] = 0.0759301442672741*xB10*dwdx19 - 0.0759301442672741*xB3*dwdx19;
    xBdot[11] = 0.0759301442672741*xB11*dwdx20 - 0.0759301442672741*xB4*dwdx20;
    xBdot[12] = -0.0759301442672741*xB10*dwdx22 - xB12*(0.0 - 0.227790432801822*dwdx21 - 0.0759301442672741*dwdx22 - 0.151860288534548*dwdx23 - 0.151860288534548*dwdx24 - 0.0759301442672741*dwdx25 - 0.0759301442672741*dwdx26 - 0.0759301442672741*dwdx27) - xB13*(0.0 - 0.0759301442672741*dwdx23 - 0.151860288534548*dwdx25 - 0.0759301442672741*dwdx27) - 0.0759301442672741*xB14*dwdx24 - 0.0759301442672741*xB15*dwdx23 - 0.0759301442672741*xB16*dwdx25 - 0.0759301442672741*xB18*dwdx27 - 0.0759301442672741*xB20*dwdx26 - xB5*(0.0 - 0.0759301442672741*dwdx24 - 0.151860288534548*dwdx26 - 0.0759301442672741*dwdx27) - 0.0759301442672741*xB7*dwdx21;
    xBdot[13] = -0.0759301442672741*xB11*dwdx29 - xB12*(0.0 - 0.151860288534548*dwdx30 - 0.0759301442672741*dwdx31 - 0.0759301442672741*dwdx34) - xB13*(0.0 - 0.227790432801822*dwdx28 - 0.0759301442672741*dwdx29 - 0.0759301442672741*dwdx30 - 0.151860288534548*dwdx31 - 0.151860288534548*dwdx32 - 0.0759301442672741*dwdx33 - 0.0759301442672741*dwdx34) - 0.0759301442672741*xB15*dwdx30 - 0.0759301442672741*xB16*dwdx31 - 0.0759301442672741*xB17*dwdx32 - 0.0759301442672741*xB18*dwdx34 - 0.0759301442672741*xB19*dwdx33 - xB5*(0.0 - 0.0759301442672741*dwdx32 - 0.151860288534548*dwdx33 - 0.0759301442672741*dwdx34) - 0.0759301442672741*xB8*dwdx28;
    xBdot[14] = -0.0759301442672741*xB10*dwdx35 - 0.0759301442672741*xB12*dwdx35 + 0.0759301442672741*xB14*dwdx35 - 0.0759301442672741*xB5*dwdx35;
    xBdot[15] = -0.0759301442672741*xB10*dwdx36 - 0.0759301442672741*xB11*dwdx37 - xB12*(0.0 + 0.0759301442672741*dwdx36 + 0.151860288534548*dwdx37) - 0.0759301442672741*xB13*dwdx36 - (-0.0759301442672741*dwdx36 - 0.0759301442672741*dwdx37)*xB15;
    xBdot[16] = -0.0759301442672741*xB10*dwdx38 - 0.0759301442672741*xB11*dwdx39 - 0.0759301442672741*xB12*dwdx39 - xB13*(0.0 + 0.151860288534548*dwdx38 + 0.0759301442672741*dwdx39) - 0.0759301442672741*xB21*dwdx40 - 0.0759301442672741*xB22*dwdx41 - 0.0759301442672741*xB23*dwdx42 - 0.0759301442672741*xB24*dwdx43 - 0.0759301442672741*xB25*dwdx44 - 0.0759301442672741*xB26*dwdx45 - 0.0759301442672741*xB27*dwdx46 - 0.0759301442672741*xB28*dwdx47 - 0.0759301442672741*xB29*dwdx48 - 0.0759301442672741*xB30*dwdx49 - 0.0759301442672741*xB31*dwdx50 - 0.0759301442672741*xB32*dwdx51 - (-0.0759301442672741*dwdx38 - 0.0759301442672741*dwdx39)*xB16;
    xBdot[17] = -0.0759301442672741*xB11*dwdx52 - 0.0759301442672741*xB13*dwdx52 + 0.0759301442672741*xB17*dwdx52 - 0.0759301442672741*xB5*dwdx52;
    xBdot[18] = -0.0759301442672741*xB10*dwdx53 - 0.0759301442672741*xB11*dwdx54 - 0.0759301442672741*xB12*dwdx54 - 0.0759301442672741*xB13*dwdx53 - 0.0759301442672741*xB21*dwdx55 - 0.0759301442672741*xB22*dwdx56 - 0.0759301442672741*xB23*dwdx57 - 0.0759301442672741*xB24*dwdx58 - 0.0759301442672741*xB25*dwdx59 - 0.0759301442672741*xB26*dwdx60 - 0.0759301442672741*xB27*dwdx61 - 0.0759301442672741*xB28*dwdx62 - 0.0759301442672741*xB29*dwdx63 - 0.0759301442672741*xB30*dwdx64 - 0.0759301442672741*xB31*dwdx65 - 0.0759301442672741*xB32*dwdx66 - (-0.0759301442672741*dwdx53 - 0.0759301442672741*dwdx54)*xB18 - (0.0759301442672741*dwdx53 + 0.0759301442672741*dwdx54)*xB5;
    xBdot[19] = -0.0759301442672741*xB11*dwdx67 + 0.0759301442672741*xB19*dwdx67 - 0.151860288534548*xB5*dwdx67;
    xBdot[20] = -0.0759301442672741*xB10*dwdx68 + 0.0759301442672741*xB20*dwdx68 - 0.0759301442672741*xB21*dwdx69 - 0.0759301442672741*xB22*dwdx70 - 0.0759301442672741*xB23*dwdx71 - 0.0759301442672741*xB24*dwdx72 - 0.0759301442672741*xB25*dwdx73 - 0.0759301442672741*xB26*dwdx74 - 0.0759301442672741*xB27*dwdx75 - 0.0759301442672741*xB28*dwdx76 - 0.0759301442672741*xB29*dwdx77 - 0.0759301442672741*xB30*dwdx78 - 0.0759301442672741*xB31*dwdx79 - 0.0759301442672741*xB32*dwdx80 - 0.151860288534548*xB5*dwdx68;
    xBdot[21] = 0.0759301442672741*xB21*dwdx81;
    xBdot[22] = 0.0759301442672741*xB22*dwdx82;
    xBdot[23] = 0.0759301442672741*xB23*dwdx83;
    xBdot[24] = 0.0759301442672741*xB24*dwdx84;
    xBdot[25] = 0.0759301442672741*xB25*dwdx85;
    xBdot[26] = 0.0759301442672741*xB26*dwdx86;
    xBdot[27] = 0.0759301442672741*xB27*dwdx87;
    xBdot[28] = 0.0759301442672741*xB28*dwdx88;
    xBdot[29] = 0.0759301442672741*xB29*dwdx89;
    xBdot[30] = 0.0759301442672741*xB30*dwdx90;
    xBdot[31] = 0.0759301442672741*xB31*dwdx91;
    xBdot[32] = 0.0759301442672741*xB32*dwdx92;
}