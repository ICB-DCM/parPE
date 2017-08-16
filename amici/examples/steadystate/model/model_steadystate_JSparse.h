#ifndef _am_model_steadystate_JSparse_h
#define _am_model_steadystate_JSparse_h

#include <sundials/sundials_direct.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_sparse.h>
#include <sundials/sundials_types.h>

class UserData;
class ReturnData;
class TempData;
class ExpData;

int JSparse_model_steadystate(realtype t, N_Vector x, N_Vector xdot, SlsMat J,
                              void *user_data, N_Vector tmp1, N_Vector tmp2,
                              N_Vector tmp3);

#endif /* _am_model_steadystate_JSparse_h */
