#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
using amici::realtype;
#include <cmath> 


#include "species.h"
#include "sensitivity.h"
#include "parameter.h"
#include "fixed_parameter.h"
#include "flux.h"
#include "dxdotdp.h"
#include "dwdx.h"
#include "JSparse.h"

void sxdot_lucarelli_12(realtype *sxdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const realtype *sx, const realtype *w, const realtype *dwdx, const realtype *J, const realtype *dxdotdp){
    sxdot[0] = dxdotdp0 + J0*sx0 + J3*sx1;
    sxdot[1] = dxdotdp1 + J1*sx0 + J4*sx1;
    sxdot[2] = dxdotdp2 + J2*sx0 + J5*sx1 + J6*sx2;
    sxdot[3] = dxdotdp3 + J11*sx3 + J17*sx5 + J28*sx6 + J39*sx10 + J7*sx2;
    sxdot[4] = dxdotdp4 + J15*sx4 + J41*sx11 + J8*sx2;
    sxdot[5] = dxdotdp5 + J12*sx3 + J18*sx5 + J29*sx6 + J37*sx9 + J43*sx12 + J53*sx13 + J63*sx14 + J89*sx17 + J93*sx18 + sx19*J111 + sx20*J114;
    sxdot[6] = dxdotdp6 + J13*sx3 + J19*sx5 + J30*sx6;
    sxdot[7] = dxdotdp7 + J31*sx7 + J44*sx12;
    sxdot[8] = dxdotdp8 + J34*sx8 + J54*sx13;
    sxdot[9] = dxdotdp9 + J20*sx5 + J38*sx9;
    sxdot[10] = dxdotdp10 + J32*sx7 + J40*sx10 + J45*sx12 + J64*sx14 + J67*sx15 + J72*sx16 + J94*sx18 + sx20*J115;
    sxdot[11] = dxdotdp11 + J35*sx8 + J42*sx11 + J55*sx13 + J68*sx15 + J73*sx16 + J90*sx17 + J95*sx18 + sx19*J112;
    sxdot[12] = dxdotdp12 + J14*sx3 + J21*sx5 + J33*sx7 + J46*sx12 + J56*sx13 + J65*sx14 + J69*sx15 + J74*sx16 + J9*sx2 + J96*sx18;
    sxdot[13] = dxdotdp13 + J10*sx2 + J16*sx4 + J22*sx5 + J36*sx8 + J47*sx12 + J57*sx13 + J70*sx15 + J75*sx16 + J91*sx17 + J97*sx18;
    sxdot[14] = dxdotdp14 + J23*sx5 + J48*sx12 + J66*sx14;
    sxdot[15] = dxdotdp15 + J49*sx12 + J58*sx13 + J71*sx15;
    sxdot[16] = dxdotdp16 + J50*sx12 + J59*sx13 + J76*sx16;
    sxdot[17] = dxdotdp17 + J24*sx5 + J60*sx13 + J92*sx17;
    sxdot[18] = dxdotdp18 + J25*sx5 + J51*sx12 + J61*sx13 + J98*sx18;
    sxdot[19] = dxdotdp19 + J26*sx5 + J62*sx13 + sx19*J113;
    sxdot[20] = dxdotdp20 + J27*sx5 + J52*sx12 + sx20*J116;
    sxdot[21] = dxdotdp21 + J77*sx16 + J99*sx18 + sx20*J117 + sx21*J129;
    sxdot[22] = dxdotdp22 + J78*sx16 + sx18*J100 + sx20*J118 + sx22*J130;
    sxdot[23] = dxdotdp23 + J79*sx16 + sx18*J101 + sx20*J119 + sx23*J131;
    sxdot[24] = dxdotdp24 + J80*sx16 + sx18*J102 + sx20*J120 + sx24*J132;
    sxdot[25] = dxdotdp25 + J81*sx16 + sx18*J103 + sx20*J121 + sx25*J133;
    sxdot[26] = dxdotdp26 + J82*sx16 + sx18*J104 + sx20*J122 + sx26*J134;
    sxdot[27] = dxdotdp27 + J83*sx16 + sx18*J105 + sx20*J123 + sx27*J135;
    sxdot[28] = dxdotdp28 + J84*sx16 + sx18*J106 + sx20*J124 + sx28*J136;
    sxdot[29] = dxdotdp29 + J85*sx16 + sx18*J107 + sx20*J125 + sx29*J137;
    sxdot[30] = dxdotdp30 + J86*sx16 + sx18*J108 + sx20*J126 + sx30*J138;
    sxdot[31] = dxdotdp31 + J87*sx16 + sx18*J109 + sx20*J127 + sx31*J139;
    sxdot[32] = dxdotdp32 + J88*sx16 + sx18*J110 + sx20*J128 + sx32*J140;
}