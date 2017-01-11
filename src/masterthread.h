#ifndef MASTER_THREAD_H
#define MASTER_THREAD_H

#include <pthread.h>

void startParameterEstimation();
void startObjectiveFunctionGradientCheck();

void *newMultiStartOptimization(void *multiStartIndexVP);
void *newLocalOptimization(void *idVP);

#endif
