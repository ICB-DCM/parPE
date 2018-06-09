#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
using amici::realtype;
#include <cmath> 

#include <sundials/sundials_sparse.h>

#include "species.h"
#include "parameter.h"
#include "fixed_parameter.h"
#include "flux.h"
#include "dwdx.h"

void JSparse_model_steadystate_scaled(SlsMat JSparse, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    JSparse->indexvals[0] = 0;
    JSparse->indexvals[1] = 1;
    JSparse->indexvals[2] = 2;
    JSparse->indexvals[3] = 0;
    JSparse->indexvals[4] = 1;
    JSparse->indexvals[5] = 2;
    JSparse->indexvals[6] = 0;
    JSparse->indexvals[7] = 1;
    JSparse->indexvals[8] = 2;
    JSparse->indexptrs[0] = 0;
    JSparse->indexptrs[1] = 3;
    JSparse->indexptrs[2] = 6;
    JSparse->indexptrs[3] = 9;
    JSparse->data[0] = 0.0 - 2.0*dwdx0 - 1.0*dwdx1;
    JSparse->data[1] = 0.0 + 1.0*dwdx0 - 1.0*dwdx1;
    JSparse->data[2] = 1.0*dwdx1;
    JSparse->data[3] = 0.0 - 1.0*dwdx2 + 2.0*dwdx3;
    JSparse->data[4] = 0.0 - 1.0*dwdx2 - 1.0*dwdx3;
    JSparse->data[5] = 1.0*dwdx2;
    JSparse->data[6] = 1.0*dwdx4;
    JSparse->data[7] = 1.0*dwdx4;
    JSparse->data[8] = -1.0*dwdx4 - 1.0*dwdx5;
}