
#include "model_steadystate_w.h"
#include <include/amici.h>
#include <include/amici_model.h>
#include <include/edata.h>
#include <include/rdata.h>
#include <include/symbolic_functions.h>
#include <include/tdata.h>
#include <include/udata.h>
#include <string.h>

int dJrzdz_model_steadystate(realtype t, int ie, N_Vector x, TempData *tdata,
                             const ExpData *edata, ReturnData *rdata) {
    int status = 0;
    Model *model = (Model *)tdata->model;
    UserData *udata = (UserData *)tdata->udata;
    realtype *x_tmp = N_VGetArrayPointer(x);
    memset(tdata->dJrzdz, 0,
           sizeof(realtype) * model->nz * model->nztrue * model->nJ);
    status = w_model_steadystate(t, x, NULL, tdata);
    return (status);
}
