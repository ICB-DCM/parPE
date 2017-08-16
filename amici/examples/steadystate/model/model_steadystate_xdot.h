#ifndef _am_model_steadystate_xdot_h
#define _am_model_steadystate_xdot_h

#include <sundials/sundials_direct.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_sparse.h>
#include <sundials/sundials_types.h>

class UserData;
class ReturnData;
class TempData;
class ExpData;

int xdot_model_steadystate(realtype t, N_Vector x, N_Vector xdot,
                           void *user_data);

#endif /* _am_model_steadystate_xdot_h */
