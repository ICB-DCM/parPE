#ifndef _am_model_steadystate_JBandB_h
#define _am_model_steadystate_JBandB_h

#include <sundials/sundials_direct.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_sparse.h>
#include <sundials/sundials_types.h>

class UserData;
class ReturnData;
class TempData;
class ExpData;

int JBandB_model_steadystate(long int NeqBdot, long int mupper, long int mlower,
                             realtype t, N_Vector x, N_Vector xB,
                             N_Vector xBdot, DlsMat JB, void *user_data,
                             N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B);

#endif /* _am_model_steadystate_JBandB_h */
