#ifndef MASTER_THREAD_H
#define MASTER_THREAD_H

#include <pthread.h>
#include <misc.h>

void startParameterEstimation(optimizerEnum optimizer);
void startObjectiveFunctionGradientCheck();

/**
 * @brief newMultiStartOptimization
 * @param multiStartIndexVP
 * @return always returns 0
 */
void *newMultiStartOptimization(void *pOptions);

/**
 * @brief newLocalOptimization
 * @param idVP
 * @return Pointer to int indicating status. 0: success, != 0: failure
 */
void *newLocalOptimization(void *pOptions);

#endif
