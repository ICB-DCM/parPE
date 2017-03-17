#ifndef CPP_MISC_H
#define CPP_MISC_H

#include <stdlib.h>

//void printMatlabArray(const double *buffer, int len);

void createDirectoryIfNotExists(char *dirName);

void strFormatCurrentLocaltime(char *buffer, size_t bufferSize, const char *format);

void shuffle(int *array, size_t numElements);

void runInParallelAndWaitForFinish(void *(*function)(void *), void **args, int numArgs);

void printBacktrace(int depth);

#endif
