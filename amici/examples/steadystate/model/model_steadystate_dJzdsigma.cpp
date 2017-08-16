
#include "model_steadystate_w.h"
#include <include/amici.h>
#include <include/edata.h>
#include <include/rdata.h>
#include <include/symbolic_functions.h>
#include <include/tdata.h>
#include <include/udata.h>
#include <string.h>

int dJzdsigma_model_steadystate(realtype t, int ie, N_Vector x, void *user_data,
                                TempData *tdata, const ExpData *edata,
                                ReturnData *rdata) {
    int status = 0;
    UserData *udata = (UserData *)user_data;
    realtype *x_tmp = N_VGetArrayPointer(x);
    memset(tdata->dJzdsigma, 0,
           sizeof(realtype) * udata->nztrue * udata->nz * udata->nJ);
    status = w_model_steadystate(t, x, NULL, user_data);
    return (status);
}
