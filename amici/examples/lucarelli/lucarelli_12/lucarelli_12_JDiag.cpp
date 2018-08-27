#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
using amici::realtype;
#include <cmath> 


#include "species.h"
#include "parameter.h"
#include "fixed_parameter.h"
#include "flux.h"
#include "dwdx.h"

void JDiag_lucarelli_12(realtype *JDiag, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    JDiag[0] = -0.0759301442672741*dwdx0;
    JDiag[1] = -0.0759301442672741*dwdx1;
    JDiag[2] = -0.0759301442672741*dwdx2;
    JDiag[3] = 0.0 - 0.0759301442672741*dwdx5 - 0.0759301442672741*dwdx6;
    JDiag[4] = -0.0759301442672741*dwdx7;
    JDiag[5] = 0.0 - 0.0759301442672741*dwdx10 - 0.0759301442672741*dwdx11 - 0.151860288534548*dwdx12 - 0.151860288534548*dwdx13 - 0.0759301442672741*dwdx14 - 0.151860288534548*dwdx8 - 0.227790432801822*dwdx9;
    JDiag[6] = -0.0759301442672741*dwdx15;
    JDiag[7] = -0.0759301442672741*dwdx16;
    JDiag[8] = -0.0759301442672741*dwdx17;
    JDiag[9] = -0.0759301442672741*dwdx18;
    JDiag[10] = -0.0759301442672741*dwdx19;
    JDiag[11] = -0.0759301442672741*dwdx20;
    JDiag[12] = 0.0 - 0.227790432801822*dwdx21 - 0.0759301442672741*dwdx22 - 0.151860288534548*dwdx23 - 0.151860288534548*dwdx24 - 0.0759301442672741*dwdx25 - 0.0759301442672741*dwdx26 - 0.0759301442672741*dwdx27;
    JDiag[13] = 0.0 - 0.227790432801822*dwdx28 - 0.0759301442672741*dwdx29 - 0.0759301442672741*dwdx30 - 0.151860288534548*dwdx31 - 0.151860288534548*dwdx32 - 0.0759301442672741*dwdx33 - 0.0759301442672741*dwdx34;
    JDiag[14] = -0.0759301442672741*dwdx35;
    JDiag[15] = -0.0759301442672741*dwdx36 - 0.0759301442672741*dwdx37;
    JDiag[16] = -0.0759301442672741*dwdx38 - 0.0759301442672741*dwdx39;
    JDiag[17] = -0.0759301442672741*dwdx52;
    JDiag[18] = -0.0759301442672741*dwdx53 - 0.0759301442672741*dwdx54;
    JDiag[19] = -0.0759301442672741*dwdx67;
    JDiag[20] = -0.0759301442672741*dwdx68;
    JDiag[21] = -0.0759301442672741*dwdx81;
    JDiag[22] = -0.0759301442672741*dwdx82;
    JDiag[23] = -0.0759301442672741*dwdx83;
    JDiag[24] = -0.0759301442672741*dwdx84;
    JDiag[25] = -0.0759301442672741*dwdx85;
    JDiag[26] = -0.0759301442672741*dwdx86;
    JDiag[27] = -0.0759301442672741*dwdx87;
    JDiag[28] = -0.0759301442672741*dwdx88;
    JDiag[29] = -0.0759301442672741*dwdx89;
    JDiag[30] = -0.0759301442672741*dwdx90;
    JDiag[31] = -0.0759301442672741*dwdx91;
    JDiag[32] = -0.0759301442672741*dwdx92;
}