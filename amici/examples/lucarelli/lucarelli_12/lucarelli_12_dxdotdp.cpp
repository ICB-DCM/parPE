#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
using amici::realtype;
#include <cmath> 


#include "species.h"
#include "parameter.h"
#include "fixed_parameter.h"
#include "flux.h"
#include "dwdp.h"

void dxdotdp_lucarelli_12(realtype *dxdotdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const realtype *w, const realtype *dwdp){
    switch(ip) {
        case 0:
            dxdotdp[0] = -0.0759301442672741*dwdp0;
            dxdotdp[1] = -0.0759301442672741*dwdp0;
            dxdotdp[2] = 0.0759301442672741*dwdp0;
            break;
        case 1:
            dxdotdp[3] = 0.0759301442672741*dwdp1;
            dxdotdp[4] = 0.0759301442672741*dwdp2;
            dxdotdp[10] = -0.0759301442672741*dwdp1;
            dxdotdp[11] = -0.0759301442672741*dwdp2;
            break;
        case 2:
            dxdotdp[5] = 0.0 + 0.0759301442672741*dwdp12 + 0.151860288534548*dwdp13 + 0.151860288534548*dwdp14 + 0.0759301442672741*dwdp15 + 0.0759301442672741*dwdp16 + 0.0759301442672741*dwdp9;
            dxdotdp[7] = -0.0759301442672741*dwdp3;
            dxdotdp[8] = -0.0759301442672741*dwdp4;
            dxdotdp[10] = 0.0 + 0.0759301442672741*dwdp10 + 0.0759301442672741*dwdp13 + 0.0759301442672741*dwdp15 + 0.0759301442672741*dwdp3 + 0.0759301442672741*dwdp5 + 0.0759301442672741*dwdp7 + 0.0759301442672741*dwdp9;
            dxdotdp[11] = 0.0 + 0.0759301442672741*dwdp11 + 0.0759301442672741*dwdp12 + 0.0759301442672741*dwdp14 + 0.0759301442672741*dwdp16 + 0.0759301442672741*dwdp4 + 0.0759301442672741*dwdp6 + 0.0759301442672741*dwdp8;
            dxdotdp[12] = 0.0 + 0.0759301442672741*dwdp11 + 0.0759301442672741*dwdp16 + 0.151860288534548*dwdp3 - 0.0759301442672741*dwdp5 + 0.0759301442672741*dwdp7 + 0.151860288534548*dwdp8 + 0.0759301442672741*dwdp9;
            dxdotdp[13] = 0.0 + 0.151860288534548*dwdp10 + 0.0759301442672741*dwdp11 + 0.0759301442672741*dwdp12 + 0.0759301442672741*dwdp15 + 0.151860288534548*dwdp4 - 0.0759301442672741*dwdp6 + 0.0759301442672741*dwdp7;
            dxdotdp[14] = -0.0759301442672741*dwdp9;
            dxdotdp[15] = -0.0759301442672741*dwdp7 - 0.0759301442672741*dwdp8;
            dxdotdp[16] = -0.0759301442672741*dwdp10 - 0.0759301442672741*dwdp11;
            dxdotdp[17] = -0.0759301442672741*dwdp12;
            dxdotdp[18] = -0.0759301442672741*dwdp15 - 0.0759301442672741*dwdp16;
            dxdotdp[19] = -0.0759301442672741*dwdp14;
            dxdotdp[20] = -0.0759301442672741*dwdp13;
            break;
        case 3:
            dxdotdp[3] = -0.0759301442672741*dwdp17;
            dxdotdp[4] = -0.0759301442672741*dwdp18;
            dxdotdp[12] = 0.0759301442672741*dwdp17;
            dxdotdp[13] = 0.0759301442672741*dwdp18;
            break;
        case 4:
            dxdotdp[21] = 0.0759301442672741*dwdp19;
            break;
        case 5:
            dxdotdp[21] = 0.0759301442672741*dwdp20;
            break;
        case 6:
            dxdotdp[21] = 0.0759301442672741*dwdp21;
            break;
        case 7:
            dxdotdp[21] = 0.0759301442672741*dwdp22;
            break;
        case 8:
            dxdotdp[21] = 0.0759301442672741*dwdp23;
            break;
        case 9:
            dxdotdp[21] = 0.0759301442672741*dwdp24;
            break;
        case 10:
            dxdotdp[21] = 0.0759301442672741*dwdp25 - 0.0759301442672741*dwdp26;
            break;
        case 11:
            dxdotdp[22] = 0.0759301442672741*dwdp27;
            break;
        case 12:
            dxdotdp[22] = 0.0759301442672741*dwdp28;
            break;
        case 13:
            dxdotdp[22] = 0.0759301442672741*dwdp29;
            break;
        case 14:
            dxdotdp[22] = 0.0759301442672741*dwdp30;
            break;
        case 15:
            dxdotdp[22] = 0.0759301442672741*dwdp31;
            break;
        case 16:
            dxdotdp[22] = 0.0759301442672741*dwdp32;
            break;
        case 17:
            dxdotdp[22] = 0.0759301442672741*dwdp33 - 0.0759301442672741*dwdp34;
            break;
        case 18:
            dxdotdp[23] = 0.0759301442672741*dwdp35;
            break;
        case 19:
            dxdotdp[23] = 0.0759301442672741*dwdp36;
            break;
        case 20:
            dxdotdp[23] = 0.0759301442672741*dwdp37;
            break;
        case 21:
            dxdotdp[23] = 0.0759301442672741*dwdp38;
            break;
        case 22:
            dxdotdp[23] = 0.0759301442672741*dwdp39;
            break;
        case 23:
            dxdotdp[23] = 0.0759301442672741*dwdp40;
            break;
        case 24:
            dxdotdp[23] = 0.0759301442672741*dwdp41 - 0.0759301442672741*dwdp42;
            break;
        case 25:
            dxdotdp[24] = 0.0759301442672741*dwdp43;
            break;
        case 26:
            dxdotdp[24] = 0.0759301442672741*dwdp44;
            break;
        case 27:
            dxdotdp[24] = 0.0759301442672741*dwdp45;
            break;
        case 28:
            dxdotdp[24] = 0.0759301442672741*dwdp46;
            break;
        case 29:
            dxdotdp[24] = 0.0759301442672741*dwdp47;
            break;
        case 30:
            dxdotdp[24] = 0.0759301442672741*dwdp48;
            break;
        case 31:
            dxdotdp[24] = 0.0759301442672741*dwdp49 - 0.0759301442672741*dwdp50;
            break;
        case 32:
            dxdotdp[25] = 0.0759301442672741*dwdp51;
            break;
        case 33:
            dxdotdp[25] = 0.0759301442672741*dwdp52;
            break;
        case 34:
            dxdotdp[25] = 0.0759301442672741*dwdp53;
            break;
        case 35:
            dxdotdp[25] = 0.0759301442672741*dwdp54;
            break;
        case 36:
            dxdotdp[25] = 0.0759301442672741*dwdp55;
            break;
        case 37:
            dxdotdp[25] = 0.0759301442672741*dwdp56;
            break;
        case 38:
            dxdotdp[25] = 0.0759301442672741*dwdp57 - 0.0759301442672741*dwdp58;
            break;
        case 39:
            dxdotdp[26] = 0.0759301442672741*dwdp59;
            break;
        case 40:
            dxdotdp[26] = 0.0759301442672741*dwdp60;
            break;
        case 41:
            dxdotdp[26] = 0.0759301442672741*dwdp61;
            break;
        case 42:
            dxdotdp[26] = 0.0759301442672741*dwdp62;
            break;
        case 43:
            dxdotdp[26] = 0.0759301442672741*dwdp63;
            break;
        case 44:
            dxdotdp[26] = 0.0759301442672741*dwdp64;
            break;
        case 45:
            dxdotdp[26] = 0.0759301442672741*dwdp65 - 0.0759301442672741*dwdp66;
            break;
        case 46:
            dxdotdp[27] = 0.0759301442672741*dwdp67;
            break;
        case 47:
            dxdotdp[27] = 0.0759301442672741*dwdp68;
            break;
        case 48:
            dxdotdp[27] = 0.0759301442672741*dwdp69;
            break;
        case 49:
            dxdotdp[27] = 0.0759301442672741*dwdp70;
            break;
        case 50:
            dxdotdp[27] = 0.0759301442672741*dwdp71;
            break;
        case 51:
            dxdotdp[27] = 0.0759301442672741*dwdp72;
            break;
        case 52:
            dxdotdp[27] = 0.0759301442672741*dwdp73 - 0.0759301442672741*dwdp74;
            break;
        case 53:
            dxdotdp[28] = 0.0759301442672741*dwdp75;
            break;
        case 54:
            dxdotdp[28] = 0.0759301442672741*dwdp76;
            break;
        case 55:
            dxdotdp[28] = 0.0759301442672741*dwdp77;
            break;
        case 56:
            dxdotdp[28] = 0.0759301442672741*dwdp78;
            break;
        case 57:
            dxdotdp[28] = 0.0759301442672741*dwdp79;
            break;
        case 58:
            dxdotdp[28] = 0.0759301442672741*dwdp80;
            break;
        case 59:
            dxdotdp[28] = 0.0759301442672741*dwdp81 - 0.0759301442672741*dwdp82;
            break;
        case 60:
            dxdotdp[29] = 0.0759301442672741*dwdp83;
            break;
        case 61:
            dxdotdp[29] = 0.0759301442672741*dwdp84;
            break;
        case 62:
            dxdotdp[29] = 0.0759301442672741*dwdp85;
            break;
        case 63:
            dxdotdp[29] = 0.0759301442672741*dwdp86;
            break;
        case 64:
            dxdotdp[29] = 0.0759301442672741*dwdp87;
            break;
        case 65:
            dxdotdp[29] = 0.0759301442672741*dwdp88;
            break;
        case 66:
            dxdotdp[29] = 0.0759301442672741*dwdp89 - 0.0759301442672741*dwdp90;
            break;
        case 67:
            dxdotdp[30] = 0.0759301442672741*dwdp91;
            break;
        case 68:
            dxdotdp[30] = 0.0759301442672741*dwdp92;
            break;
        case 69:
            dxdotdp[30] = 0.0759301442672741*dwdp93;
            break;
        case 70:
            dxdotdp[30] = 0.0759301442672741*dwdp94;
            break;
        case 71:
            dxdotdp[30] = 0.0759301442672741*dwdp95;
            break;
        case 72:
            dxdotdp[30] = 0.0759301442672741*dwdp96;
            break;
        case 73:
            dxdotdp[30] = 0.0759301442672741*dwdp97 - 0.0759301442672741*dwdp98;
            break;
        case 74:
            dxdotdp[31] = 0.0759301442672741*dwdp99;
            break;
        case 75:
            dxdotdp[31] = 0.0759301442672741*dwdp100;
            break;
        case 76:
            dxdotdp[31] = 0.0759301442672741*dwdp101;
            break;
        case 77:
            dxdotdp[31] = 0.0759301442672741*dwdp102;
            break;
        case 78:
            dxdotdp[31] = 0.0759301442672741*dwdp103;
            break;
        case 79:
            dxdotdp[31] = 0.0759301442672741*dwdp104;
            break;
        case 80:
            dxdotdp[31] = 0.0759301442672741*dwdp105 - 0.0759301442672741*dwdp106;
            break;
        case 81:
            dxdotdp[32] = 0.0759301442672741*dwdp107;
            break;
        case 82:
            dxdotdp[32] = 0.0759301442672741*dwdp108;
            break;
        case 83:
            dxdotdp[32] = 0.0759301442672741*dwdp109;
            break;
        case 84:
            dxdotdp[32] = 0.0759301442672741*dwdp110;
            break;
        case 85:
            dxdotdp[32] = 0.0759301442672741*dwdp111;
            break;
        case 86:
            dxdotdp[32] = 0.0759301442672741*dwdp112;
            break;
        case 87:
            dxdotdp[32] = 0.0759301442672741*dwdp113 - 0.0759301442672741*dwdp114;
            break;
        case 88:
            break;
        case 89:
            dxdotdp[12] = -0.151860288534548*dwdp115;
            dxdotdp[13] = -0.0759301442672741*dwdp115;
            dxdotdp[15] = 0.0759301442672741*dwdp115;
            break;
        case 90:
            dxdotdp[5] = -0.0759301442672741*dwdp116;
            dxdotdp[12] = -0.151860288534548*dwdp116;
            dxdotdp[14] = 0.0759301442672741*dwdp116;
            break;
        case 91:
            dxdotdp[12] = -0.0759301442672741*dwdp117;
            dxdotdp[13] = -0.151860288534548*dwdp117;
            dxdotdp[16] = 0.0759301442672741*dwdp117;
            break;
        case 92:
            dxdotdp[5] = -0.0759301442672741*dwdp118;
            dxdotdp[12] = -0.0759301442672741*dwdp118;
            dxdotdp[13] = -0.0759301442672741*dwdp118;
            dxdotdp[18] = 0.0759301442672741*dwdp118;
            break;
        case 93:
            dxdotdp[5] = -0.151860288534548*dwdp119;
            dxdotdp[12] = -0.0759301442672741*dwdp119;
            dxdotdp[20] = 0.0759301442672741*dwdp119;
            break;
        case 94:
            dxdotdp[5] = -0.0759301442672741*dwdp120;
            dxdotdp[13] = -0.151860288534548*dwdp120;
            dxdotdp[17] = 0.0759301442672741*dwdp120;
            break;
        case 95:
            dxdotdp[5] = -0.151860288534548*dwdp121;
            dxdotdp[13] = -0.0759301442672741*dwdp121;
            dxdotdp[19] = 0.0759301442672741*dwdp121;
            break;
        case 96:
            dxdotdp[3] = -0.0759301442672741*dwdp122;
            dxdotdp[5] = -0.151860288534548*dwdp122;
            dxdotdp[6] = 0.0759301442672741*dwdp122;
            break;
        case 97:
            dxdotdp[3] = 0.0759301442672741*dwdp123;
            dxdotdp[5] = 0.0 + 0.151860288534548*dwdp123 + 0.227790432801822*dwdp124;
            dxdotdp[6] = -0.0759301442672741*dwdp123;
            dxdotdp[9] = -0.0759301442672741*dwdp124;
            break;
        case 98:
            dxdotdp[7] = 0.0759301442672741*dwdp125;
            dxdotdp[12] = -0.227790432801822*dwdp125;
            break;
        case 99:
            dxdotdp[8] = 0.0759301442672741*dwdp126;
            dxdotdp[13] = -0.227790432801822*dwdp126;
            break;
        case 100:
            dxdotdp[5] = -0.227790432801822*dwdp127;
            dxdotdp[9] = 0.0759301442672741*dwdp127;
            break;
        case 101:
            dxdotdp[2] = -0.0759301442672741*dwdp128;
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