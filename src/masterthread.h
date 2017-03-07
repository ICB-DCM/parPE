#ifndef MASTER_THREAD_H
#define MASTER_THREAD_H

#include <pthread.h>

void startParameterEstimation();
void startObjectiveFunctionGradientCheck();

/**
 * @brief newMultiStartOptimization
 * @param multiStartIndexVP
 * @return always returns 0
 */
void *newMultiStartOptimization(void *multiStartIndexVP);

/**
 * @brief newLocalOptimization
 * @param idVP
 * @return Pointer to int indicating status. 0: success, != 0: failure
 */
void *newLocalOptimization(void *idVP);

#endif
