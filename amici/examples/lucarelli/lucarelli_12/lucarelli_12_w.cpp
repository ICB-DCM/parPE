#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
using amici::realtype;
#include <cmath> 


#include "species.h"
#include "parameter.h"
#include "fixed_parameter.h"

void w_lucarelli_12(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h){
    w[0] = 13.17*Rec*TGFb*Rec_act;
    w[1] = 13.17*TGFb_pRec*pRec_degind;
    w[2] = 13.17*S2*pow(S4, 2)*k_on_u;
    w[3] = 13.17*S2_S4_S4*kdiss_SS;
    w[4] = 13.17*pow(ppS2, 3)*khomo2;
    w[5] = 39.51*ppS2_ppS2_ppS2*S_dephosphos;
    w[6] = 13.17*pow(ppS3, 3)*khomo3;
    w[7] = 39.51*S_dephosphos*ppS3_ppS3_ppS3;
    w[8] = 13.17*pow(S4, 3)*khomo4;
    w[9] = 13.17*S4_S4_S4*kdiss_SS;
    w[10] = 13.17*S2*TGFb_pRec*S_phos;
    w[11] = 13.17*ppS2*S_dephosphos;
    w[12] = 13.17*pS2*S_dephos;
    w[13] = 13.17*S3*TGFb_pRec*S_phos;
    w[14] = 13.17*ppS3*S_dephosphos;
    w[15] = 13.17*pS3*S_dephos;
    w[16] = 13.17*pow(ppS2, 2)*ppS3*k_223;
    w[17] = 26.34*ppS2_ppS2_ppS3*S_dephosphos;
    w[18] = 13.17*ppS2_ppS2_ppS3*S_dephosphos;
    w[19] = 13.17*S4*pow(ppS2, 2)*k_224;
    w[20] = 26.34*S_dephosphos*ppS2_ppS2_S4;
    w[21] = 13.17*ppS2*pow(ppS3, 2)*k_233;
    w[22] = 13.17*ppS2_ppS3_ppS3*S_dephosphos;
    w[23] = 26.34*ppS2_ppS3_ppS3*S_dephosphos;
    w[24] = 13.17*S4*pow(ppS3, 2)*k_334;
    w[25] = 26.34*S_dephosphos*ppS3_ppS3_S4;
    w[26] = 13.17*pow(S4, 2)*ppS2*k_244;
    w[27] = 13.17*S_dephosphos*ppS2_S4_S4;
    w[28] = 13.17*pow(S4, 2)*ppS3*k_344;
    w[29] = 13.17*S_dephosphos*ppS3_S4_S4;
    w[30] = 13.17*S4*ppS2*ppS3*k_234;
    w[31] = 13.17*S_dephosphos*ppS2_ppS3_S4;
    w[32] = 13.17*S_dephosphos*ppS2_ppS3_S4;
    w[33] = 13.17*(geneA_turn + ppS2_S4_S4*geneA_act2 + ppS2_ppS3_S4*geneA_act1 + ppS2_ppS3_ppS3*geneA_act3)/(1 + ppS2_S4_S4*geneA_inh2 + ppS2_ppS3_S4*geneA_inh1 + ppS2_ppS3_ppS3*geneA_inh3);
    w[34] = 13.17*geneA*geneA_turn;
    w[35] = 13.17*(geneB_turn + ppS2_S4_S4*geneB_act2 + ppS2_ppS3_S4*geneB_act1 + ppS2_ppS3_ppS3*geneB_act3)/(1 + ppS2_S4_S4*geneB_inh2 + ppS2_ppS3_S4*geneB_inh1 + ppS2_ppS3_ppS3*geneB_inh3);
    w[36] = 13.17*geneB*geneB_turn;
    w[37] = 13.17*(geneC_turn + ppS2_S4_S4*geneC_act2 + ppS2_ppS3_S4*geneC_act1 + ppS2_ppS3_ppS3*geneC_act3)/(1 + ppS2_S4_S4*geneC_inh2 + ppS2_ppS3_S4*geneC_inh1 + ppS2_ppS3_ppS3*geneC_inh3);
    w[38] = 13.17*geneC*geneC_turn;
    w[39] = 13.17*(geneD_turn + ppS2_S4_S4*geneD_act2 + ppS2_ppS3_S4*geneD_act1 + ppS2_ppS3_ppS3*geneD_act3)/(1 + ppS2_S4_S4*geneD_inh2 + ppS2_ppS3_S4*geneD_inh1 + ppS2_ppS3_ppS3*geneD_inh3);
    w[40] = 13.17*geneD*geneD_turn;
    w[41] = 13.17*(geneE_turn + ppS2_S4_S4*geneE_act2 + ppS2_ppS3_S4*geneE_act1 + ppS2_ppS3_ppS3*geneE_act3)/(1 + ppS2_S4_S4*geneE_inh2 + ppS2_ppS3_S4*geneE_inh1 + ppS2_ppS3_ppS3*geneE_inh3);
    w[42] = 13.17*geneE*geneE_turn;
    w[43] = 13.17*(geneF_turn + ppS2_S4_S4*geneF_act2 + ppS2_ppS3_S4*geneF_act1 + ppS2_ppS3_ppS3*geneF_act3)/(1 + ppS2_S4_S4*geneF_inh2 + ppS2_ppS3_S4*geneF_inh1 + ppS2_ppS3_ppS3*geneF_inh3);
    w[44] = 13.17*geneF*geneF_turn;
    w[45] = 13.17*(geneG_turn + ppS2_S4_S4*geneG_act2 + ppS2_ppS3_S4*geneG_act1 + ppS2_ppS3_ppS3*geneG_act3)/(1 + ppS2_S4_S4*geneG_inh2 + ppS2_ppS3_S4*geneG_inh1 + ppS2_ppS3_ppS3*geneG_inh3);
    w[46] = 13.17*geneG*geneG_turn;
    w[47] = 13.17*(geneH_turn + ppS2_S4_S4*geneH_act2 + ppS2_ppS3_S4*geneH_act1 + ppS2_ppS3_ppS3*geneH_act3)/(1 + ppS2_S4_S4*geneH_inh2 + ppS2_ppS3_S4*geneH_inh1 + ppS2_ppS3_ppS3*geneH_inh3);
    w[48] = 13.17*geneH*geneH_turn;
    w[49] = 13.17*(geneI_turn + ppS2_S4_S4*geneI_act2 + ppS2_ppS3_S4*geneI_act1 + ppS2_ppS3_ppS3*geneI_act3)/(1 + ppS2_S4_S4*geneI_inh2 + ppS2_ppS3_S4*geneI_inh1 + ppS2_ppS3_ppS3*geneI_inh3);
    w[50] = 13.17*geneI*geneI_turn;
    w[51] = 13.17*(geneJ_turn + ppS2_S4_S4*geneJ_act2 + ppS2_ppS3_S4*geneJ_act1 + ppS2_ppS3_ppS3*geneJ_act3)/(1 + ppS2_S4_S4*geneJ_inh2 + ppS2_ppS3_S4*geneJ_inh1 + ppS2_ppS3_ppS3*geneJ_inh3);
    w[52] = 13.17*geneJ*geneJ_turn;
    w[53] = 13.17*(geneK_turn + ppS2_S4_S4*geneK_act2 + ppS2_ppS3_S4*geneK_act1 + ppS2_ppS3_ppS3*geneK_act3)/(1 + ppS2_S4_S4*geneK_inh2 + ppS2_ppS3_S4*geneK_inh1 + ppS2_ppS3_ppS3*geneK_inh3);
    w[54] = 13.17*geneK*geneK_turn;
    w[55] = 13.17*(geneL_turn + ppS2_S4_S4*geneL_act2 + ppS2_ppS3_S4*geneL_act1 + ppS2_ppS3_ppS3*geneL_act3)/(1 + ppS2_S4_S4*geneL_inh2 + ppS2_ppS3_S4*geneL_inh1 + ppS2_ppS3_ppS3*geneL_inh3);
    w[56] = 13.17*geneL*geneL_turn;
}