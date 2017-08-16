
#include "model_steadystate_w.h"
#include <include/amici.h>
#include <include/edata.h>
#include <include/rdata.h>
#include <include/symbolic_functions.h>
#include <include/tdata.h>
#include <include/udata.h>
#include <string.h>

int dJydsigma_model_steadystate(realtype t, int it, N_Vector x, void *user_data,
                                TempData *tdata, const ExpData *edata,
                                ReturnData *rdata) {
    int status = 0;
    UserData *udata = (UserData *)user_data;
    realtype *x_tmp = N_VGetArrayPointer(x);
    memset(tdata->dJydsigma, 0,
           sizeof(realtype) * udata->nytrue * udata->ny * udata->nJ);
    status = w_model_steadystate(t, x, NULL, user_data);
    int iy;
    if (!amiIsNaN(edata->my[0 * udata->nt + it])) {
        iy = 0;
        tdata->dJydsigma[iy + (0 + 0 * 1) * udata->nytrue] =
            1.0 / (tdata->sigmay[0] * tdata->sigmay[0] * tdata->sigmay[0]) *
                pow(edata->my[it + udata->nt * 0] -
                        rdata->y[it + udata->nt * 0],
                    2.0) *
                -1.0 +
            1.0 / tdata->sigmay[0];
    }
    if (!amiIsNaN(edata->my[1 * udata->nt + it])) {
        iy = 1;
        tdata->dJydsigma[iy + (0 + 1 * 1) * udata->nytrue] =
            1.0 / (tdata->sigmay[1] * tdata->sigmay[1] * tdata->sigmay[1]) *
                pow(edata->my[it + udata->nt * 1] -
                        rdata->y[it + udata->nt * 1],
                    2.0) *
                -1.0 +
            1.0 / tdata->sigmay[1];
    }
    if (!amiIsNaN(edata->my[2 * udata->nt + it])) {
        iy = 2;
        tdata->dJydsigma[iy + (0 + 2 * 1) * udata->nytrue] =
            1.0 / (tdata->sigmay[2] * tdata->sigmay[2] * tdata->sigmay[2]) *
                pow(edata->my[it + udata->nt * 2] -
                        rdata->y[it + udata->nt * 2],
                    2.0) *
                -1.0 +
            1.0 / tdata->sigmay[2];
    }
    return (status);
}
