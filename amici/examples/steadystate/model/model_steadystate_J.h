#ifndef _am_model_steadystate_J_h
#define _am_model_steadystate_J_h

#include <sundials/sundials_direct.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_sparse.h>
#include <sundials/sundials_types.h>

class UserData;
class ReturnData;
class TempData;
class ExpData;

int J_model_steadystate(long int N, realtype t, N_Vector x, N_Vector xdot,
                        DlsMat J, void *user_data, N_Vector tmp1, N_Vector tmp2,
                        N_Vector tmp3);

#endif /* _am_model_steadystate_J_h */
