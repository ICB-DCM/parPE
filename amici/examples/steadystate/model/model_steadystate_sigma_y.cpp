
#include "model_steadystate_w.h"
#include <include/amici.h>
#include <include/symbolic_functions.h>
#include <include/tdata.h>
#include <include/udata.h>
#include <string.h>

int sigma_y_model_steadystate(realtype t, void *user_data, TempData *tdata) {
    int status = 0;
    UserData *udata = (UserData *)user_data;
    memset(tdata->sigmay, 0, sizeof(realtype) * 3);
    tdata->sigmay[0] = 1.0;
    tdata->sigmay[1] = 1.0;
    tdata->sigmay[2] = 1.0;
    return (status);
}
