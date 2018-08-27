#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
using amici::realtype;
#include <cmath> 


#include "species.h"
#include "parameter.h"
#include "fixed_parameter.h"
#include "flux.h"

void dwdx_lucarelli_12(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    dwdx[0] = 13.17*Rec*Rec_act;
    dwdx[1] = 13.17*TGFb*Rec_act;
    dwdx[2] = 13.17*pRec_degind;
    dwdx[3] = 13.17*S2*S_phos;
    dwdx[4] = 13.17*S3*S_phos;
    dwdx[5] = 13.17*pow(S4, 2)*k_on_u;
    dwdx[6] = 13.17*TGFb_pRec*S_phos;
    dwdx[7] = 13.17*TGFb_pRec*S_phos;
    dwdx[8] = 26.34*S2*S4*k_on_u;
    dwdx[9] = 39.51*pow(S4, 2)*khomo4;
    dwdx[10] = 13.17*pow(ppS2, 2)*k_224;
    dwdx[11] = 13.17*pow(ppS3, 2)*k_334;
    dwdx[12] = 26.34*S4*ppS2*k_244;
    dwdx[13] = 26.34*S4*ppS3*k_344;
    dwdx[14] = 13.17*ppS2*ppS3*k_234;
    dwdx[15] = 13.17*kdiss_SS;
    dwdx[16] = 39.51*S_dephosphos;
    dwdx[17] = 39.51*S_dephosphos;
    dwdx[18] = 13.17*kdiss_SS;
    dwdx[19] = 13.17*S_dephos;
    dwdx[20] = 13.17*S_dephos;
    dwdx[21] = 39.51*pow(ppS2, 2)*khomo2;
    dwdx[22] = 13.17*S_dephosphos;
    dwdx[23] = 26.34*ppS2*ppS3*k_223;
    dwdx[24] = 26.34*S4*ppS2*k_224;
    dwdx[25] = 13.17*pow(ppS3, 2)*k_233;
    dwdx[26] = 13.17*pow(S4, 2)*k_244;
    dwdx[27] = 13.17*S4*ppS3*k_234;
    dwdx[28] = 39.51*pow(ppS3, 2)*khomo3;
    dwdx[29] = 13.17*S_dephosphos;
    dwdx[30] = 13.17*pow(ppS2, 2)*k_223;
    dwdx[31] = 26.34*ppS2*ppS3*k_233;
    dwdx[32] = 26.34*S4*ppS3*k_334;
    dwdx[33] = 13.17*pow(S4, 2)*k_344;
    dwdx[34] = 13.17*S4*ppS2*k_234;
    dwdx[35] = 26.34*S_dephosphos;
    dwdx[36] = 26.34*S_dephosphos;
    dwdx[37] = 13.17*S_dephosphos;
    dwdx[38] = 13.17*S_dephosphos;
    dwdx[39] = 26.34*S_dephosphos;
    dwdx[40] = 13.17*geneA_act3/(1 + ppS2_S4_S4*geneA_inh2 + ppS2_ppS3_S4*geneA_inh1 + ppS2_ppS3_ppS3*geneA_inh3) - 13.17*(geneA_turn + ppS2_S4_S4*geneA_act2 + ppS2_ppS3_S4*geneA_act1 + ppS2_ppS3_ppS3*geneA_act3)*geneA_inh3/pow(1 + ppS2_S4_S4*geneA_inh2 + ppS2_ppS3_S4*geneA_inh1 + ppS2_ppS3_ppS3*geneA_inh3, 2);
    dwdx[41] = 13.17*geneB_act3/(1 + ppS2_S4_S4*geneB_inh2 + ppS2_ppS3_S4*geneB_inh1 + ppS2_ppS3_ppS3*geneB_inh3) - 13.17*(geneB_turn + ppS2_S4_S4*geneB_act2 + ppS2_ppS3_S4*geneB_act1 + ppS2_ppS3_ppS3*geneB_act3)*geneB_inh3/pow(1 + ppS2_S4_S4*geneB_inh2 + ppS2_ppS3_S4*geneB_inh1 + ppS2_ppS3_ppS3*geneB_inh3, 2);
    dwdx[42] = 13.17*geneC_act3/(1 + ppS2_S4_S4*geneC_inh2 + ppS2_ppS3_S4*geneC_inh1 + ppS2_ppS3_ppS3*geneC_inh3) - 13.17*(geneC_turn + ppS2_S4_S4*geneC_act2 + ppS2_ppS3_S4*geneC_act1 + ppS2_ppS3_ppS3*geneC_act3)*geneC_inh3/pow(1 + ppS2_S4_S4*geneC_inh2 + ppS2_ppS3_S4*geneC_inh1 + ppS2_ppS3_ppS3*geneC_inh3, 2);
    dwdx[43] = 13.17*geneD_act3/(1 + ppS2_S4_S4*geneD_inh2 + ppS2_ppS3_S4*geneD_inh1 + ppS2_ppS3_ppS3*geneD_inh3) - 13.17*(geneD_turn + ppS2_S4_S4*geneD_act2 + ppS2_ppS3_S4*geneD_act1 + ppS2_ppS3_ppS3*geneD_act3)*geneD_inh3/pow(1 + ppS2_S4_S4*geneD_inh2 + ppS2_ppS3_S4*geneD_inh1 + ppS2_ppS3_ppS3*geneD_inh3, 2);
    dwdx[44] = 13.17*geneE_act3/(1 + ppS2_S4_S4*geneE_inh2 + ppS2_ppS3_S4*geneE_inh1 + ppS2_ppS3_ppS3*geneE_inh3) - 13.17*(geneE_turn + ppS2_S4_S4*geneE_act2 + ppS2_ppS3_S4*geneE_act1 + ppS2_ppS3_ppS3*geneE_act3)*geneE_inh3/pow(1 + ppS2_S4_S4*geneE_inh2 + ppS2_ppS3_S4*geneE_inh1 + ppS2_ppS3_ppS3*geneE_inh3, 2);
    dwdx[45] = 13.17*geneF_act3/(1 + ppS2_S4_S4*geneF_inh2 + ppS2_ppS3_S4*geneF_inh1 + ppS2_ppS3_ppS3*geneF_inh3) - 13.17*(geneF_turn + ppS2_S4_S4*geneF_act2 + ppS2_ppS3_S4*geneF_act1 + ppS2_ppS3_ppS3*geneF_act3)*geneF_inh3/pow(1 + ppS2_S4_S4*geneF_inh2 + ppS2_ppS3_S4*geneF_inh1 + ppS2_ppS3_ppS3*geneF_inh3, 2);
    dwdx[46] = 13.17*geneG_act3/(1 + ppS2_S4_S4*geneG_inh2 + ppS2_ppS3_S4*geneG_inh1 + ppS2_ppS3_ppS3*geneG_inh3) - 13.17*(geneG_turn + ppS2_S4_S4*geneG_act2 + ppS2_ppS3_S4*geneG_act1 + ppS2_ppS3_ppS3*geneG_act3)*geneG_inh3/pow(1 + ppS2_S4_S4*geneG_inh2 + ppS2_ppS3_S4*geneG_inh1 + ppS2_ppS3_ppS3*geneG_inh3, 2);
    dwdx[47] = 13.17*geneH_act3/(1 + ppS2_S4_S4*geneH_inh2 + ppS2_ppS3_S4*geneH_inh1 + ppS2_ppS3_ppS3*geneH_inh3) - 13.17*(geneH_turn + ppS2_S4_S4*geneH_act2 + ppS2_ppS3_S4*geneH_act1 + ppS2_ppS3_ppS3*geneH_act3)*geneH_inh3/pow(1 + ppS2_S4_S4*geneH_inh2 + ppS2_ppS3_S4*geneH_inh1 + ppS2_ppS3_ppS3*geneH_inh3, 2);
    dwdx[48] = 13.17*geneI_act3/(1 + ppS2_S4_S4*geneI_inh2 + ppS2_ppS3_S4*geneI_inh1 + ppS2_ppS3_ppS3*geneI_inh3) - 13.17*(geneI_turn + ppS2_S4_S4*geneI_act2 + ppS2_ppS3_S4*geneI_act1 + ppS2_ppS3_ppS3*geneI_act3)*geneI_inh3/pow(1 + ppS2_S4_S4*geneI_inh2 + ppS2_ppS3_S4*geneI_inh1 + ppS2_ppS3_ppS3*geneI_inh3, 2);
    dwdx[49] = 13.17*geneJ_act3/(1 + ppS2_S4_S4*geneJ_inh2 + ppS2_ppS3_S4*geneJ_inh1 + ppS2_ppS3_ppS3*geneJ_inh3) - 13.17*(geneJ_turn + ppS2_S4_S4*geneJ_act2 + ppS2_ppS3_S4*geneJ_act1 + ppS2_ppS3_ppS3*geneJ_act3)*geneJ_inh3/pow(1 + ppS2_S4_S4*geneJ_inh2 + ppS2_ppS3_S4*geneJ_inh1 + ppS2_ppS3_ppS3*geneJ_inh3, 2);
    dwdx[50] = 13.17*geneK_act3/(1 + ppS2_S4_S4*geneK_inh2 + ppS2_ppS3_S4*geneK_inh1 + ppS2_ppS3_ppS3*geneK_inh3) - 13.17*(geneK_turn + ppS2_S4_S4*geneK_act2 + ppS2_ppS3_S4*geneK_act1 + ppS2_ppS3_ppS3*geneK_act3)*geneK_inh3/pow(1 + ppS2_S4_S4*geneK_inh2 + ppS2_ppS3_S4*geneK_inh1 + ppS2_ppS3_ppS3*geneK_inh3, 2);
    dwdx[51] = 13.17*geneL_act3/(1 + ppS2_S4_S4*geneL_inh2 + ppS2_ppS3_S4*geneL_inh1 + ppS2_ppS3_ppS3*geneL_inh3) - 13.17*(geneL_turn + ppS2_S4_S4*geneL_act2 + ppS2_ppS3_S4*geneL_act1 + ppS2_ppS3_ppS3*geneL_act3)*geneL_inh3/pow(1 + ppS2_S4_S4*geneL_inh2 + ppS2_ppS3_S4*geneL_inh1 + ppS2_ppS3_ppS3*geneL_inh3, 2);
    dwdx[52] = 26.34*S_dephosphos;
    dwdx[53] = 13.17*S_dephosphos;
    dwdx[54] = 13.17*S_dephosphos;
    dwdx[55] = 13.17*geneA_act1/(1 + ppS2_S4_S4*geneA_inh2 + ppS2_ppS3_S4*geneA_inh1 + ppS2_ppS3_ppS3*geneA_inh3) - 13.17*(geneA_turn + ppS2_S4_S4*geneA_act2 + ppS2_ppS3_S4*geneA_act1 + ppS2_ppS3_ppS3*geneA_act3)*geneA_inh1/pow(1 + ppS2_S4_S4*geneA_inh2 + ppS2_ppS3_S4*geneA_inh1 + ppS2_ppS3_ppS3*geneA_inh3, 2);
    dwdx[56] = 13.17*geneB_act1/(1 + ppS2_S4_S4*geneB_inh2 + ppS2_ppS3_S4*geneB_inh1 + ppS2_ppS3_ppS3*geneB_inh3) - 13.17*(geneB_turn + ppS2_S4_S4*geneB_act2 + ppS2_ppS3_S4*geneB_act1 + ppS2_ppS3_ppS3*geneB_act3)*geneB_inh1/pow(1 + ppS2_S4_S4*geneB_inh2 + ppS2_ppS3_S4*geneB_inh1 + ppS2_ppS3_ppS3*geneB_inh3, 2);
    dwdx[57] = 13.17*geneC_act1/(1 + ppS2_S4_S4*geneC_inh2 + ppS2_ppS3_S4*geneC_inh1 + ppS2_ppS3_ppS3*geneC_inh3) - 13.17*(geneC_turn + ppS2_S4_S4*geneC_act2 + ppS2_ppS3_S4*geneC_act1 + ppS2_ppS3_ppS3*geneC_act3)*geneC_inh1/pow(1 + ppS2_S4_S4*geneC_inh2 + ppS2_ppS3_S4*geneC_inh1 + ppS2_ppS3_ppS3*geneC_inh3, 2);
    dwdx[58] = 13.17*geneD_act1/(1 + ppS2_S4_S4*geneD_inh2 + ppS2_ppS3_S4*geneD_inh1 + ppS2_ppS3_ppS3*geneD_inh3) - 13.17*(geneD_turn + ppS2_S4_S4*geneD_act2 + ppS2_ppS3_S4*geneD_act1 + ppS2_ppS3_ppS3*geneD_act3)*geneD_inh1/pow(1 + ppS2_S4_S4*geneD_inh2 + ppS2_ppS3_S4*geneD_inh1 + ppS2_ppS3_ppS3*geneD_inh3, 2);
    dwdx[59] = 13.17*geneE_act1/(1 + ppS2_S4_S4*geneE_inh2 + ppS2_ppS3_S4*geneE_inh1 + ppS2_ppS3_ppS3*geneE_inh3) - 13.17*(geneE_turn + ppS2_S4_S4*geneE_act2 + ppS2_ppS3_S4*geneE_act1 + ppS2_ppS3_ppS3*geneE_act3)*geneE_inh1/pow(1 + ppS2_S4_S4*geneE_inh2 + ppS2_ppS3_S4*geneE_inh1 + ppS2_ppS3_ppS3*geneE_inh3, 2);
    dwdx[60] = 13.17*geneF_act1/(1 + ppS2_S4_S4*geneF_inh2 + ppS2_ppS3_S4*geneF_inh1 + ppS2_ppS3_ppS3*geneF_inh3) - 13.17*(geneF_turn + ppS2_S4_S4*geneF_act2 + ppS2_ppS3_S4*geneF_act1 + ppS2_ppS3_ppS3*geneF_act3)*geneF_inh1/pow(1 + ppS2_S4_S4*geneF_inh2 + ppS2_ppS3_S4*geneF_inh1 + ppS2_ppS3_ppS3*geneF_inh3, 2);
    dwdx[61] = 13.17*geneG_act1/(1 + ppS2_S4_S4*geneG_inh2 + ppS2_ppS3_S4*geneG_inh1 + ppS2_ppS3_ppS3*geneG_inh3) - 13.17*(geneG_turn + ppS2_S4_S4*geneG_act2 + ppS2_ppS3_S4*geneG_act1 + ppS2_ppS3_ppS3*geneG_act3)*geneG_inh1/pow(1 + ppS2_S4_S4*geneG_inh2 + ppS2_ppS3_S4*geneG_inh1 + ppS2_ppS3_ppS3*geneG_inh3, 2);
    dwdx[62] = 13.17*geneH_act1/(1 + ppS2_S4_S4*geneH_inh2 + ppS2_ppS3_S4*geneH_inh1 + ppS2_ppS3_ppS3*geneH_inh3) - 13.17*(geneH_turn + ppS2_S4_S4*geneH_act2 + ppS2_ppS3_S4*geneH_act1 + ppS2_ppS3_ppS3*geneH_act3)*geneH_inh1/pow(1 + ppS2_S4_S4*geneH_inh2 + ppS2_ppS3_S4*geneH_inh1 + ppS2_ppS3_ppS3*geneH_inh3, 2);
    dwdx[63] = 13.17*geneI_act1/(1 + ppS2_S4_S4*geneI_inh2 + ppS2_ppS3_S4*geneI_inh1 + ppS2_ppS3_ppS3*geneI_inh3) - 13.17*(geneI_turn + ppS2_S4_S4*geneI_act2 + ppS2_ppS3_S4*geneI_act1 + ppS2_ppS3_ppS3*geneI_act3)*geneI_inh1/pow(1 + ppS2_S4_S4*geneI_inh2 + ppS2_ppS3_S4*geneI_inh1 + ppS2_ppS3_ppS3*geneI_inh3, 2);
    dwdx[64] = 13.17*geneJ_act1/(1 + ppS2_S4_S4*geneJ_inh2 + ppS2_ppS3_S4*geneJ_inh1 + ppS2_ppS3_ppS3*geneJ_inh3) - 13.17*(geneJ_turn + ppS2_S4_S4*geneJ_act2 + ppS2_ppS3_S4*geneJ_act1 + ppS2_ppS3_ppS3*geneJ_act3)*geneJ_inh1/pow(1 + ppS2_S4_S4*geneJ_inh2 + ppS2_ppS3_S4*geneJ_inh1 + ppS2_ppS3_ppS3*geneJ_inh3, 2);
    dwdx[65] = 13.17*geneK_act1/(1 + ppS2_S4_S4*geneK_inh2 + ppS2_ppS3_S4*geneK_inh1 + ppS2_ppS3_ppS3*geneK_inh3) - 13.17*(geneK_turn + ppS2_S4_S4*geneK_act2 + ppS2_ppS3_S4*geneK_act1 + ppS2_ppS3_ppS3*geneK_act3)*geneK_inh1/pow(1 + ppS2_S4_S4*geneK_inh2 + ppS2_ppS3_S4*geneK_inh1 + ppS2_ppS3_ppS3*geneK_inh3, 2);
    dwdx[66] = 13.17*geneL_act1/(1 + ppS2_S4_S4*geneL_inh2 + ppS2_ppS3_S4*geneL_inh1 + ppS2_ppS3_ppS3*geneL_inh3) - 13.17*(geneL_turn + ppS2_S4_S4*geneL_act2 + ppS2_ppS3_S4*geneL_act1 + ppS2_ppS3_ppS3*geneL_act3)*geneL_inh1/pow(1 + ppS2_S4_S4*geneL_inh2 + ppS2_ppS3_S4*geneL_inh1 + ppS2_ppS3_ppS3*geneL_inh3, 2);
    dwdx[67] = 13.17*S_dephosphos;
    dwdx[68] = 13.17*S_dephosphos;
    dwdx[69] = 13.17*geneA_act2/(1 + ppS2_S4_S4*geneA_inh2 + ppS2_ppS3_S4*geneA_inh1 + ppS2_ppS3_ppS3*geneA_inh3) - 13.17*(geneA_turn + ppS2_S4_S4*geneA_act2 + ppS2_ppS3_S4*geneA_act1 + ppS2_ppS3_ppS3*geneA_act3)*geneA_inh2/pow(1 + ppS2_S4_S4*geneA_inh2 + ppS2_ppS3_S4*geneA_inh1 + ppS2_ppS3_ppS3*geneA_inh3, 2);
    dwdx[70] = 13.17*geneB_act2/(1 + ppS2_S4_S4*geneB_inh2 + ppS2_ppS3_S4*geneB_inh1 + ppS2_ppS3_ppS3*geneB_inh3) - 13.17*(geneB_turn + ppS2_S4_S4*geneB_act2 + ppS2_ppS3_S4*geneB_act1 + ppS2_ppS3_ppS3*geneB_act3)*geneB_inh2/pow(1 + ppS2_S4_S4*geneB_inh2 + ppS2_ppS3_S4*geneB_inh1 + ppS2_ppS3_ppS3*geneB_inh3, 2);
    dwdx[71] = 13.17*geneC_act2/(1 + ppS2_S4_S4*geneC_inh2 + ppS2_ppS3_S4*geneC_inh1 + ppS2_ppS3_ppS3*geneC_inh3) - 13.17*(geneC_turn + ppS2_S4_S4*geneC_act2 + ppS2_ppS3_S4*geneC_act1 + ppS2_ppS3_ppS3*geneC_act3)*geneC_inh2/pow(1 + ppS2_S4_S4*geneC_inh2 + ppS2_ppS3_S4*geneC_inh1 + ppS2_ppS3_ppS3*geneC_inh3, 2);
    dwdx[72] = 13.17*geneD_act2/(1 + ppS2_S4_S4*geneD_inh2 + ppS2_ppS3_S4*geneD_inh1 + ppS2_ppS3_ppS3*geneD_inh3) - 13.17*(geneD_turn + ppS2_S4_S4*geneD_act2 + ppS2_ppS3_S4*geneD_act1 + ppS2_ppS3_ppS3*geneD_act3)*geneD_inh2/pow(1 + ppS2_S4_S4*geneD_inh2 + ppS2_ppS3_S4*geneD_inh1 + ppS2_ppS3_ppS3*geneD_inh3, 2);
    dwdx[73] = 13.17*geneE_act2/(1 + ppS2_S4_S4*geneE_inh2 + ppS2_ppS3_S4*geneE_inh1 + ppS2_ppS3_ppS3*geneE_inh3) - 13.17*(geneE_turn + ppS2_S4_S4*geneE_act2 + ppS2_ppS3_S4*geneE_act1 + ppS2_ppS3_ppS3*geneE_act3)*geneE_inh2/pow(1 + ppS2_S4_S4*geneE_inh2 + ppS2_ppS3_S4*geneE_inh1 + ppS2_ppS3_ppS3*geneE_inh3, 2);
    dwdx[74] = 13.17*geneF_act2/(1 + ppS2_S4_S4*geneF_inh2 + ppS2_ppS3_S4*geneF_inh1 + ppS2_ppS3_ppS3*geneF_inh3) - 13.17*(geneF_turn + ppS2_S4_S4*geneF_act2 + ppS2_ppS3_S4*geneF_act1 + ppS2_ppS3_ppS3*geneF_act3)*geneF_inh2/pow(1 + ppS2_S4_S4*geneF_inh2 + ppS2_ppS3_S4*geneF_inh1 + ppS2_ppS3_ppS3*geneF_inh3, 2);
    dwdx[75] = 13.17*geneG_act2/(1 + ppS2_S4_S4*geneG_inh2 + ppS2_ppS3_S4*geneG_inh1 + ppS2_ppS3_ppS3*geneG_inh3) - 13.17*(geneG_turn + ppS2_S4_S4*geneG_act2 + ppS2_ppS3_S4*geneG_act1 + ppS2_ppS3_ppS3*geneG_act3)*geneG_inh2/pow(1 + ppS2_S4_S4*geneG_inh2 + ppS2_ppS3_S4*geneG_inh1 + ppS2_ppS3_ppS3*geneG_inh3, 2);
    dwdx[76] = 13.17*geneH_act2/(1 + ppS2_S4_S4*geneH_inh2 + ppS2_ppS3_S4*geneH_inh1 + ppS2_ppS3_ppS3*geneH_inh3) - 13.17*(geneH_turn + ppS2_S4_S4*geneH_act2 + ppS2_ppS3_S4*geneH_act1 + ppS2_ppS3_ppS3*geneH_act3)*geneH_inh2/pow(1 + ppS2_S4_S4*geneH_inh2 + ppS2_ppS3_S4*geneH_inh1 + ppS2_ppS3_ppS3*geneH_inh3, 2);
    dwdx[77] = 13.17*geneI_act2/(1 + ppS2_S4_S4*geneI_inh2 + ppS2_ppS3_S4*geneI_inh1 + ppS2_ppS3_ppS3*geneI_inh3) - 13.17*(geneI_turn + ppS2_S4_S4*geneI_act2 + ppS2_ppS3_S4*geneI_act1 + ppS2_ppS3_ppS3*geneI_act3)*geneI_inh2/pow(1 + ppS2_S4_S4*geneI_inh2 + ppS2_ppS3_S4*geneI_inh1 + ppS2_ppS3_ppS3*geneI_inh3, 2);
    dwdx[78] = 13.17*geneJ_act2/(1 + ppS2_S4_S4*geneJ_inh2 + ppS2_ppS3_S4*geneJ_inh1 + ppS2_ppS3_ppS3*geneJ_inh3) - 13.17*(geneJ_turn + ppS2_S4_S4*geneJ_act2 + ppS2_ppS3_S4*geneJ_act1 + ppS2_ppS3_ppS3*geneJ_act3)*geneJ_inh2/pow(1 + ppS2_S4_S4*geneJ_inh2 + ppS2_ppS3_S4*geneJ_inh1 + ppS2_ppS3_ppS3*geneJ_inh3, 2);
    dwdx[79] = 13.17*geneK_act2/(1 + ppS2_S4_S4*geneK_inh2 + ppS2_ppS3_S4*geneK_inh1 + ppS2_ppS3_ppS3*geneK_inh3) - 13.17*(geneK_turn + ppS2_S4_S4*geneK_act2 + ppS2_ppS3_S4*geneK_act1 + ppS2_ppS3_ppS3*geneK_act3)*geneK_inh2/pow(1 + ppS2_S4_S4*geneK_inh2 + ppS2_ppS3_S4*geneK_inh1 + ppS2_ppS3_ppS3*geneK_inh3, 2);
    dwdx[80] = 13.17*geneL_act2/(1 + ppS2_S4_S4*geneL_inh2 + ppS2_ppS3_S4*geneL_inh1 + ppS2_ppS3_ppS3*geneL_inh3) - 13.17*(geneL_turn + ppS2_S4_S4*geneL_act2 + ppS2_ppS3_S4*geneL_act1 + ppS2_ppS3_ppS3*geneL_act3)*geneL_inh2/pow(1 + ppS2_S4_S4*geneL_inh2 + ppS2_ppS3_S4*geneL_inh1 + ppS2_ppS3_ppS3*geneL_inh3, 2);
    dwdx[81] = 13.17*geneA_turn;
    dwdx[82] = 13.17*geneB_turn;
    dwdx[83] = 13.17*geneC_turn;
    dwdx[84] = 13.17*geneD_turn;
    dwdx[85] = 13.17*geneE_turn;
    dwdx[86] = 13.17*geneF_turn;
    dwdx[87] = 13.17*geneG_turn;
    dwdx[88] = 13.17*geneH_turn;
    dwdx[89] = 13.17*geneI_turn;
    dwdx[90] = 13.17*geneJ_turn;
    dwdx[91] = 13.17*geneK_turn;
    dwdx[92] = 13.17*geneL_turn;
}