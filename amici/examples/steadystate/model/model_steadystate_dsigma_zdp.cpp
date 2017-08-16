
#include "model_steadystate_w.h"
#include <include/amici.h>
#include <include/symbolic_functions.h>
#include <include/tdata.h>
#include <include/udata.h>
#include <string.h>

int dsigma_zdp_model_steadystate(realtype t, int ie, void *user_data,
                                 TempData *tdata) {
    int status = 0;
    UserData *udata = (UserData *)user_data;
    int ip;
    memset(tdata->dsigmazdp, 0, sizeof(realtype) * 0 * udata->nplist);
    for (ip = 0; ip < udata->nplist; ip++) {
        switch (udata->plist[ip]) {}
    }
    return (status);
}
