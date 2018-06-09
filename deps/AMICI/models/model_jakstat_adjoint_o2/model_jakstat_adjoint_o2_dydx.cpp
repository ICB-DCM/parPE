
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

void dydx_model_jakstat_adjoint_o2(double *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h) {
  dydx[0+1*54] = p[13]/p[4];
  dydx[0+2*54] = (p[13]*2.0)/p[4];
  dydx[1+0*54] = p[12]/p[4];
  dydx[1+1*54] = p[12]/p[4];
  dydx[1+2*54] = (p[12]*2.0)/p[4];
  dydx[3+10*54] = p[13]/p[4];
  dydx[3+11*54] = (p[13]*2.0)/p[4];
  dydx[4+9*54] = p[12]/p[4];
  dydx[4+10*54] = p[12]/p[4];
  dydx[4+11*54] = (p[12]*2.0)/p[4];
  dydx[6+19*54] = p[13]/p[4];
  dydx[6+20*54] = (p[13]*2.0)/p[4];
  dydx[7+18*54] = p[12]/p[4];
  dydx[7+19*54] = p[12]/p[4];
  dydx[7+20*54] = (p[12]*2.0)/p[4];
  dydx[9+28*54] = p[13]/p[4];
  dydx[9+29*54] = (p[13]*2.0)/p[4];
  dydx[10+27*54] = p[12]/p[4];
  dydx[10+28*54] = p[12]/p[4];
  dydx[10+29*54] = (p[12]*2.0)/p[4];
  dydx[12+37*54] = p[13]/p[4];
  dydx[12+38*54] = (p[13]*2.0)/p[4];
  dydx[13+36*54] = p[12]/p[4];
  dydx[13+37*54] = p[12]/p[4];
  dydx[13+38*54] = (p[12]*2.0)/p[4];
  dydx[15+1*54] = -1.0/(p[4]*p[4])*p[13];
  dydx[15+2*54] = 1.0/(p[4]*p[4])*p[13]*-2.0;
  dydx[15+46*54] = p[13]/p[4];
  dydx[15+47*54] = (p[13]*2.0)/p[4];
  dydx[16+0*54] = -1.0/(p[4]*p[4])*p[12];
  dydx[16+1*54] = -1.0/(p[4]*p[4])*p[12];
  dydx[16+2*54] = 1.0/(p[4]*p[4])*p[12]*-2.0;
  dydx[16+45*54] = p[12]/p[4];
  dydx[16+46*54] = p[12]/p[4];
  dydx[16+47*54] = (p[12]*2.0)/p[4];
  dydx[18+55*54] = p[13]/p[4];
  dydx[18+56*54] = (p[13]*2.0)/p[4];
  dydx[19+54*54] = p[12]/p[4];
  dydx[19+55*54] = p[12]/p[4];
  dydx[19+56*54] = (p[12]*2.0)/p[4];
  dydx[21+64*54] = p[13]/p[4];
  dydx[21+65*54] = (p[13]*2.0)/p[4];
  dydx[22+63*54] = p[12]/p[4];
  dydx[22+64*54] = p[12]/p[4];
  dydx[22+65*54] = (p[12]*2.0)/p[4];
  dydx[24+73*54] = p[13]/p[4];
  dydx[24+74*54] = (p[13]*2.0)/p[4];
  dydx[25+72*54] = p[12]/p[4];
  dydx[25+73*54] = p[12]/p[4];
  dydx[25+74*54] = (p[12]*2.0)/p[4];
  dydx[27+82*54] = p[13]/p[4];
  dydx[27+83*54] = (p[13]*2.0)/p[4];
  dydx[28+81*54] = p[12]/p[4];
  dydx[28+82*54] = p[12]/p[4];
  dydx[28+83*54] = (p[12]*2.0)/p[4];
  dydx[30+91*54] = p[13]/p[4];
  dydx[30+92*54] = (p[13]*2.0)/p[4];
  dydx[31+90*54] = p[12]/p[4];
  dydx[31+91*54] = p[12]/p[4];
  dydx[31+92*54] = (p[12]*2.0)/p[4];
  dydx[33+100*54] = p[13]/p[4];
  dydx[33+101*54] = (p[13]*2.0)/p[4];
  dydx[34+99*54] = p[12]/p[4];
  dydx[34+100*54] = p[12]/p[4];
  dydx[34+101*54] = (p[12]*2.0)/p[4];
  dydx[36+109*54] = p[13]/p[4];
  dydx[36+110*54] = (p[13]*2.0)/p[4];
  dydx[37+108*54] = p[12]/p[4];
  dydx[37+109*54] = p[12]/p[4];
  dydx[37+110*54] = (p[12]*2.0)/p[4];
  dydx[39+118*54] = p[13]/p[4];
  dydx[39+119*54] = (p[13]*2.0)/p[4];
  dydx[40+0*54] = 1.0/p[4];
  dydx[40+1*54] = 1.0/p[4];
  dydx[40+2*54] = 2.0/p[4];
  dydx[40+117*54] = p[12]/p[4];
  dydx[40+118*54] = p[12]/p[4];
  dydx[40+119*54] = (p[12]*2.0)/p[4];
  dydx[42+1*54] = 1.0/p[4];
  dydx[42+2*54] = 2.0/p[4];
  dydx[42+127*54] = p[13]/p[4];
  dydx[42+128*54] = (p[13]*2.0)/p[4];
  dydx[43+126*54] = p[12]/p[4];
  dydx[43+127*54] = p[12]/p[4];
  dydx[43+128*54] = (p[12]*2.0)/p[4];
  dydx[45+136*54] = p[13]/p[4];
  dydx[45+137*54] = (p[13]*2.0)/p[4];
  dydx[46+135*54] = p[12]/p[4];
  dydx[46+136*54] = p[12]/p[4];
  dydx[46+137*54] = (p[12]*2.0)/p[4];
  dydx[48+145*54] = p[13]/p[4];
  dydx[48+146*54] = (p[13]*2.0)/p[4];
  dydx[49+144*54] = p[12]/p[4];
  dydx[49+145*54] = p[12]/p[4];
  dydx[49+146*54] = (p[12]*2.0)/p[4];
  dydx[51+154*54] = p[13]/p[4];
  dydx[51+155*54] = (p[13]*2.0)/p[4];
  dydx[52+153*54] = p[12]/p[4];
  dydx[52+154*54] = p[12]/p[4];
  dydx[52+155*54] = (p[12]*2.0)/p[4];
}

