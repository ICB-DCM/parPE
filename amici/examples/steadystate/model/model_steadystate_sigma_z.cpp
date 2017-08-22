
#include "model_steadystate_w.h"
#include <include/amici.h>
#include <include/amici_model.h>
#include <include/symbolic_functions.h>
#include <include/tdata.h>
#include <include/udata.h>
#include <string.h>

int sigma_z_model_steadystate(realtype t, int ie, TempData *tdata) {
    int status = 0;
    Model *model = (Model *)tdata->model;
    UserData *udata = (UserData *)tdata->udata;
    memset(tdata->sigmaz, 0, sizeof(realtype) * 0);
    return (status);
}
