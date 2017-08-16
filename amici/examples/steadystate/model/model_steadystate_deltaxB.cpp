
#include "model_steadystate_w.h"
#include <include/amici.h>
#include <include/symbolic_functions.h>
#include <include/tdata.h>
#include <include/udata.h>
#include <string.h>

int deltaxB_model_steadystate(realtype t, int ie, N_Vector x, N_Vector xB,
                              N_Vector xdot, N_Vector xdot_old, void *user_data,
                              TempData *tdata) {
    int status = 0;
    UserData *udata = (UserData *)user_data;
    realtype *x_tmp = N_VGetArrayPointer(x);
    realtype *xB_tmp = N_VGetArrayPointer(xB);
    realtype *xdot_tmp = N_VGetArrayPointer(xdot);
    realtype *xdot_old_tmp = N_VGetArrayPointer(xdot_old);
    memset(tdata->deltaxB, 0, sizeof(realtype) * 3);
    status = w_model_steadystate(t, x, NULL, user_data);
    return (status);
}
