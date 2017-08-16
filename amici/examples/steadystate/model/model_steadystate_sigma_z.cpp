
#include "model_steadystate_w.h"
#include <include/amici.h>
#include <include/symbolic_functions.h>
#include <include/tdata.h>
#include <include/udata.h>
#include <string.h>

int sigma_z_model_steadystate(realtype t, int ie, void *user_data,
                              TempData *tdata) {
    int status = 0;
    UserData *udata = (UserData *)user_data;
    memset(tdata->sigmaz, 0, sizeof(realtype) * 0);
    return (status);
}
