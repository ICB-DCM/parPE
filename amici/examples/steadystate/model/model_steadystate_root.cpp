
#include "model_steadystate_w.h"
#include <include/amici.h>
#include <include/symbolic_functions.h>
#include <include/udata.h>
#include <string.h>

int root_model_steadystate(realtype t, N_Vector x, realtype *root,
                           void *user_data) {
    int status = 0;
    UserData *udata = (UserData *)user_data;
    realtype *x_tmp = N_VGetArrayPointer(x);
    status = w_model_steadystate(t, x, NULL, user_data);
    return (status);
}
