#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
using amici::realtype;
#include <cmath> 

#include <sundials/sundials_sparse.h>

#include "species.h"
#include "parameter.h"
#include "fixed_parameter.h"
#include "adjoint.h"
#include "flux.h"
#include "dwdx.h"

void JSparseB_model_steadystate_scaled(SlsMat JSparseB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx){
    JSparseB->indexvals[0] = 0;
    JSparseB->indexvals[1] = 1;
    JSparseB->indexvals[2] = 2;
    JSparseB->indexvals[3] = 0;
    JSparseB->indexvals[4] = 1;
    JSparseB->indexvals[5] = 2;
    JSparseB->indexvals[6] = 0;
    JSparseB->indexvals[7] = 1;
    JSparseB->indexvals[8] = 2;
    JSparseB->indexptrs[0] = 0;
    JSparseB->indexptrs[1] = 3;
    JSparseB->indexptrs[2] = 6;
    JSparseB->indexptrs[3] = 9;
    JSparseB->data[0] = 0.0 - 2.0*dwdx0 - 1.0*dwdx1;
    JSparseB->data[1] = 0.0 - 1.0*dwdx2 + 2.0*dwdx3;
    JSparseB->data[2] = 1.0*dwdx4;
    JSparseB->data[3] = 0.0 + 1.0*dwdx0 - 1.0*dwdx1;
    JSparseB->data[4] = 0.0 - 1.0*dwdx2 - 1.0*dwdx3;
    JSparseB->data[5] = 1.0*dwdx4;
    JSparseB->data[6] = 1.0*dwdx1;
    JSparseB->data[7] = 1.0*dwdx2;
    JSparseB->data[8] = -1.0*dwdx4 - 1.0*dwdx5;
}