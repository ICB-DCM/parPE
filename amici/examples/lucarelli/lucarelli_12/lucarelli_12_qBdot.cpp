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

void qBdot_lucarelli_12(realtype *qBdot, const int ip, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdp){
    switch(ip) {
        case 0:
            qBdot[0] = 0.0759301442672741*xB0*dwdp0 + 0.0759301442672741*xB1*dwdp0 - 0.0759301442672741*xB2*dwdp0;
            break;
        case 1:
            qBdot[0] = 0.0759301442672741*xB10*dwdp1 + 0.0759301442672741*xB11*dwdp2 - 0.0759301442672741*xB3*dwdp1 - 0.0759301442672741*xB4*dwdp2;
            break;
        case 2:
            qBdot[0] = -xB10*(0.0 + 0.0759301442672741*dwdp10 + 0.0759301442672741*dwdp13 + 0.0759301442672741*dwdp15 + 0.0759301442672741*dwdp3 + 0.0759301442672741*dwdp5 + 0.0759301442672741*dwdp7 + 0.0759301442672741*dwdp9) - xB11*(0.0 + 0.0759301442672741*dwdp11 + 0.0759301442672741*dwdp12 + 0.0759301442672741*dwdp14 + 0.0759301442672741*dwdp16 + 0.0759301442672741*dwdp4 + 0.0759301442672741*dwdp6 + 0.0759301442672741*dwdp8) - xB12*(0.0 + 0.0759301442672741*dwdp11 + 0.0759301442672741*dwdp16 + 0.151860288534548*dwdp3 - 0.0759301442672741*dwdp5 + 0.0759301442672741*dwdp7 + 0.151860288534548*dwdp8 + 0.0759301442672741*dwdp9) - xB13*(0.0 + 0.151860288534548*dwdp10 + 0.0759301442672741*dwdp11 + 0.0759301442672741*dwdp12 + 0.0759301442672741*dwdp15 + 0.151860288534548*dwdp4 - 0.0759301442672741*dwdp6 + 0.0759301442672741*dwdp7) + 0.0759301442672741*xB14*dwdp9 + 0.0759301442672741*xB17*dwdp12 + 0.0759301442672741*xB19*dwdp14 + 0.0759301442672741*xB20*dwdp13 - xB5*(0.0 + 0.0759301442672741*dwdp12 + 0.151860288534548*dwdp13 + 0.151860288534548*dwdp14 + 0.0759301442672741*dwdp15 + 0.0759301442672741*dwdp16 + 0.0759301442672741*dwdp9) + 0.0759301442672741*xB7*dwdp3 + 0.0759301442672741*xB8*dwdp4 - (-0.0759301442672741*dwdp10 - 0.0759301442672741*dwdp11)*xB16 - (-0.0759301442672741*dwdp15 - 0.0759301442672741*dwdp16)*xB18 - (-0.0759301442672741*dwdp7 - 0.0759301442672741*dwdp8)*xB15;
            break;
        case 3:
            qBdot[0] = -0.0759301442672741*xB12*dwdp17 - 0.0759301442672741*xB13*dwdp18 + 0.0759301442672741*xB3*dwdp17 + 0.0759301442672741*xB4*dwdp18;
            break;
        case 4:
            qBdot[0] = -0.0759301442672741*xB21*dwdp19;
            break;
        case 5:
            qBdot[0] = -0.0759301442672741*xB21*dwdp20;
            break;
        case 6:
            qBdot[0] = -0.0759301442672741*xB21*dwdp21;
            break;
        case 7:
            qBdot[0] = -0.0759301442672741*xB21*dwdp22;
            break;
        case 8:
            qBdot[0] = -0.0759301442672741*xB21*dwdp23;
            break;
        case 9:
            qBdot[0] = -0.0759301442672741*xB21*dwdp24;
            break;
        case 10:
            qBdot[0] = -xB21*(0.0759301442672741*dwdp25 - 0.0759301442672741*dwdp26);
            break;
        case 11:
            qBdot[0] = -0.0759301442672741*xB22*dwdp27;
            break;
        case 12:
            qBdot[0] = -0.0759301442672741*xB22*dwdp28;
            break;
        case 13:
            qBdot[0] = -0.0759301442672741*xB22*dwdp29;
            break;
        case 14:
            qBdot[0] = -0.0759301442672741*xB22*dwdp30;
            break;
        case 15:
            qBdot[0] = -0.0759301442672741*xB22*dwdp31;
            break;
        case 16:
            qBdot[0] = -0.0759301442672741*xB22*dwdp32;
            break;
        case 17:
            qBdot[0] = -xB22*(0.0759301442672741*dwdp33 - 0.0759301442672741*dwdp34);
            break;
        case 18:
            qBdot[0] = -0.0759301442672741*xB23*dwdp35;
            break;
        case 19:
            qBdot[0] = -0.0759301442672741*xB23*dwdp36;
            break;
        case 20:
            qBdot[0] = -0.0759301442672741*xB23*dwdp37;
            break;
        case 21:
            qBdot[0] = -0.0759301442672741*xB23*dwdp38;
            break;
        case 22:
            qBdot[0] = -0.0759301442672741*xB23*dwdp39;
            break;
        case 23:
            qBdot[0] = -0.0759301442672741*xB23*dwdp40;
            break;
        case 24:
            qBdot[0] = -xB23*(0.0759301442672741*dwdp41 - 0.0759301442672741*dwdp42);
            break;
        case 25:
            qBdot[0] = -0.0759301442672741*xB24*dwdp43;
            break;
        case 26:
            qBdot[0] = -0.0759301442672741*xB24*dwdp44;
            break;
        case 27:
            qBdot[0] = -0.0759301442672741*xB24*dwdp45;
            break;
        case 28:
            qBdot[0] = -0.0759301442672741*xB24*dwdp46;
            break;
        case 29:
            qBdot[0] = -0.0759301442672741*xB24*dwdp47;
            break;
        case 30:
            qBdot[0] = -0.0759301442672741*xB24*dwdp48;
            break;
        case 31:
            qBdot[0] = -xB24*(0.0759301442672741*dwdp49 - 0.0759301442672741*dwdp50);
            break;
        case 32:
            qBdot[0] = -0.0759301442672741*xB25*dwdp51;
            break;
        case 33:
            qBdot[0] = -0.0759301442672741*xB25*dwdp52;
            break;
        case 34:
            qBdot[0] = -0.0759301442672741*xB25*dwdp53;
            break;
        case 35:
            qBdot[0] = -0.0759301442672741*xB25*dwdp54;
            break;
        case 36:
            qBdot[0] = -0.0759301442672741*xB25*dwdp55;
            break;
        case 37:
            qBdot[0] = -0.0759301442672741*xB25*dwdp56;
            break;
        case 38:
            qBdot[0] = -xB25*(0.0759301442672741*dwdp57 - 0.0759301442672741*dwdp58);
            break;
        case 39:
            qBdot[0] = -0.0759301442672741*xB26*dwdp59;
            break;
        case 40:
            qBdot[0] = -0.0759301442672741*xB26*dwdp60;
            break;
        case 41:
            qBdot[0] = -0.0759301442672741*xB26*dwdp61;
            break;
        case 42:
            qBdot[0] = -0.0759301442672741*xB26*dwdp62;
            break;
        case 43:
            qBdot[0] = -0.0759301442672741*xB26*dwdp63;
            break;
        case 44:
            qBdot[0] = -0.0759301442672741*xB26*dwdp64;
            break;
        case 45:
            qBdot[0] = -xB26*(0.0759301442672741*dwdp65 - 0.0759301442672741*dwdp66);
            break;
        case 46:
            qBdot[0] = -0.0759301442672741*xB27*dwdp67;
            break;
        case 47:
            qBdot[0] = -0.0759301442672741*xB27*dwdp68;
            break;
        case 48:
            qBdot[0] = -0.0759301442672741*xB27*dwdp69;
            break;
        case 49:
            qBdot[0] = -0.0759301442672741*xB27*dwdp70;
            break;
        case 50:
            qBdot[0] = -0.0759301442672741*xB27*dwdp71;
            break;
        case 51:
            qBdot[0] = -0.0759301442672741*xB27*dwdp72;
            break;
        case 52:
            qBdot[0] = -xB27*(0.0759301442672741*dwdp73 - 0.0759301442672741*dwdp74);
            break;
        case 53:
            qBdot[0] = -0.0759301442672741*xB28*dwdp75;
            break;
        case 54:
            qBdot[0] = -0.0759301442672741*xB28*dwdp76;
            break;
        case 55:
            qBdot[0] = -0.0759301442672741*xB28*dwdp77;
            break;
        case 56:
            qBdot[0] = -0.0759301442672741*xB28*dwdp78;
            break;
        case 57:
            qBdot[0] = -0.0759301442672741*xB28*dwdp79;
            break;
        case 58:
            qBdot[0] = -0.0759301442672741*xB28*dwdp80;
            break;
        case 59:
            qBdot[0] = -xB28*(0.0759301442672741*dwdp81 - 0.0759301442672741*dwdp82);
            break;
        case 60:
            qBdot[0] = -0.0759301442672741*xB29*dwdp83;
            break;
        case 61:
            qBdot[0] = -0.0759301442672741*xB29*dwdp84;
            break;
        case 62:
            qBdot[0] = -0.0759301442672741*xB29*dwdp85;
            break;
        case 63:
            qBdot[0] = -0.0759301442672741*xB29*dwdp86;
            break;
        case 64:
            qBdot[0] = -0.0759301442672741*xB29*dwdp87;
            break;
        case 65:
            qBdot[0] = -0.0759301442672741*xB29*dwdp88;
            break;
        case 66:
            qBdot[0] = -xB29*(0.0759301442672741*dwdp89 - 0.0759301442672741*dwdp90);
            break;
        case 67:
            qBdot[0] = -0.0759301442672741*xB30*dwdp91;
            break;
        case 68:
            qBdot[0] = -0.0759301442672741*xB30*dwdp92;
            break;
        case 69:
            qBdot[0] = -0.0759301442672741*xB30*dwdp93;
            break;
        case 70:
            qBdot[0] = -0.0759301442672741*xB30*dwdp94;
            break;
        case 71:
            qBdot[0] = -0.0759301442672741*xB30*dwdp95;
            break;
        case 72:
            qBdot[0] = -0.0759301442672741*xB30*dwdp96;
            break;
        case 73:
            qBdot[0] = -xB30*(0.0759301442672741*dwdp97 - 0.0759301442672741*dwdp98);
            break;
        case 74:
            qBdot[0] = -0.0759301442672741*xB31*dwdp99;
            break;
        case 75:
            qBdot[0] = -0.0759301442672741*xB31*dwdp100;
            break;
        case 76:
            qBdot[0] = -0.0759301442672741*xB31*dwdp101;
            break;
        case 77:
            qBdot[0] = -0.0759301442672741*xB31*dwdp102;
            break;
        case 78:
            qBdot[0] = -0.0759301442672741*xB31*dwdp103;
            break;
        case 79:
            qBdot[0] = -0.0759301442672741*xB31*dwdp104;
            break;
        case 80:
            qBdot[0] = -xB31*(0.0759301442672741*dwdp105 - 0.0759301442672741*dwdp106);
            break;
        case 81:
            qBdot[0] = -0.0759301442672741*xB32*dwdp107;
            break;
        case 82:
            qBdot[0] = -0.0759301442672741*xB32*dwdp108;
            break;
        case 83:
            qBdot[0] = -0.0759301442672741*xB32*dwdp109;
            break;
        case 84:
            qBdot[0] = -0.0759301442672741*xB32*dwdp110;
            break;
        case 85:
            qBdot[0] = -0.0759301442672741*xB32*dwdp111;
            break;
        case 86:
            qBdot[0] = -0.0759301442672741*xB32*dwdp112;
            break;
        case 87:
            qBdot[0] = -xB32*(0.0759301442672741*dwdp113 - 0.0759301442672741*dwdp114);
            break;
        case 88:
            break;
        case 89:
            qBdot[0] = 0.151860288534548*xB12*dwdp115 + 0.0759301442672741*xB13*dwdp115 - 0.0759301442672741*xB15*dwdp115;
            break;
        case 90:
            qBdot[0] = 0.151860288534548*xB12*dwdp116 - 0.0759301442672741*xB14*dwdp116 + 0.0759301442672741*xB5*dwdp116;
            break;
        case 91:
            qBdot[0] = 0.0759301442672741*xB12*dwdp117 + 0.151860288534548*xB13*dwdp117 - 0.0759301442672741*xB16*dwdp117;
            break;
        case 92:
            qBdot[0] = 0.0759301442672741*xB12*dwdp118 + 0.0759301442672741*xB13*dwdp118 - 0.0759301442672741*xB18*dwdp118 + 0.0759301442672741*xB5*dwdp118;
            break;
        case 93:
            qBdot[0] = 0.0759301442672741*xB12*dwdp119 - 0.0759301442672741*xB20*dwdp119 + 0.151860288534548*xB5*dwdp119;
            break;
        case 94:
            qBdot[0] = 0.151860288534548*xB13*dwdp120 - 0.0759301442672741*xB17*dwdp120 + 0.0759301442672741*xB5*dwdp120;
            break;
        case 95:
            qBdot[0] = 0.0759301442672741*xB13*dwdp121 - 0.0759301442672741*xB19*dwdp121 + 0.151860288534548*xB5*dwdp121;
            break;
        case 96:
            qBdot[0] = 0.0759301442672741*xB3*dwdp122 + 0.151860288534548*xB5*dwdp122 - 0.0759301442672741*xB6*dwdp122;
            break;
        case 97:
            qBdot[0] = -0.0759301442672741*xB3*dwdp123 - xB5*(0.0 + 0.151860288534548*dwdp123 + 0.227790432801822*dwdp124) + 0.0759301442672741*xB6*dwdp123 + 0.0759301442672741*xB9*dwdp124;
            break;
        case 98:
            qBdot[0] = 0.227790432801822*xB12*dwdp125 - 0.0759301442672741*xB7*dwdp125;
            break;
        case 99:
            qBdot[0] = 0.227790432801822*xB13*dwdp126 - 0.0759301442672741*xB8*dwdp126;
            break;
        case 100:
            qBdot[0] = 0.227790432801822*xB5*dwdp127 - 0.0759301442672741*xB9*dwdp127;
            break;
        case 101:
            qBdot[0] = 0.0759301442672741*xB2*dwdp128;
            break;
        case 102:
            break;
        case 103:
            break;
        case 104:
            break;
        case 105:
            break;
        case 106:
            break;
        case 107:
            break;
        case 108:
            break;
        case 109:
            break;
        case 110:
            break;
        case 111:
            break;
        case 112:
            break;
        case 113:
            break;
}
}