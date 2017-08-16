
#include "model_steadystate_w.h"
#include <include/amici.h>
#include <include/symbolic_functions.h>
#include <include/udata.h>
#include <string.h>

int dwdp_model_steadystate(realtype t, N_Vector x, N_Vector dx,
                           void *user_data) {
    int status = 0;
    UserData *udata = (UserData *)user_data;
    realtype *x_tmp = N_VGetArrayPointer(x);
    memset(udata->dwdp, 0, sizeof(realtype) * 1);
    status = w_model_steadystate(t, x, NULL, user_data);
    udata->dwdp[0] = x_tmp[2];
    return (status);
}
