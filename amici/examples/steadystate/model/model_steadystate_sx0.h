#ifndef _am_model_steadystate_sx0_h
#define _am_model_steadystate_sx0_h

#include <sundials/sundials_direct.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_sparse.h>
#include <sundials/sundials_types.h>

class UserData;
class ReturnData;
class TempData;
class ExpData;

int sx0_model_steadystate(N_Vector *sx0, N_Vector x, N_Vector dx,
                          void *user_data);

#endif /* _am_model_steadystate_sx0_h */
