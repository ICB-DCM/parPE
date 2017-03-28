#include "misc.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <time.h>
#include <alloca.h>
#include <unistd.h>
#include <math.h>
#include <sys/stat.h>
#include "logging.h"
#include <execinfo.h>

//void printMatlabArray(const double *buffer, int len)
//{
//    printf("[");
//    printfArray(buffer, len - 1, "%e, ");
//    printf("%e]\n", buffer[len - 1]);
//    fflush(stdout);
//}


void createDirectoryIfNotExists(char *dirName)
{
    struct stat st = {0};

    if (stat(dirName, &st) == -1) {
        mkdir(dirName, 0700);
    }
}

void strFormatCurrentLocaltime(char *buffer, size_t bufferSize, const char *format) {
    time_t timer;
    time(&timer);

    struct tm* tm_info;
    tm_info = localtime(&timer);

    strftime(buffer, bufferSize, format, tm_info);
}

void shuffle(int *array, size_t numElements)
{
    size_t i;
    for (i = 0; i < numElements - 1; ++i) {
        size_t j = numElements * rand() / RAND_MAX;
        int tmp = array[j];
        array[j] = array[i];
        array[i] = tmp;
    }
}


#include <pthread.h>

void runInParallelAndWaitForFinish(void *(*function)(void *), void **args, int numArgs) {
    // create threads
    pthread_attr_t threadAttr;
    pthread_attr_init(&threadAttr);
    pthread_attr_setdetachstate(&threadAttr, PTHREAD_CREATE_JOINABLE);

    pthread_t *threads = alloca(numArgs * sizeof(pthread_t));

    for(int i = 0; i < numArgs; ++i) {
        pthread_create(&threads[i], &threadAttr, function, args[i]);
    }
    pthread_attr_destroy(&threadAttr);

    // wait for finish
    for(int i = 0; i < numArgs; ++i) {
        pthread_join(threads[i], NULL);
        logmessage(LOGLVL_DEBUG, "Thread i %d finished", i);
    }
    logmessage(LOGLVL_DEBUG, "All k threads finished.");
}

void printBacktrace(int depth) {
    void *array[depth];
    size_t size;
    size = backtrace(array, depth);
    backtrace_symbols_fd(array, size, STDERR_FILENO);
}

double randDouble(double min, double max) {
    return min + rand() / (double) RAND_MAX * (max - min);
}

void fillArrayRandomDoubleIndividualInterval(const double *min, const double *max, int length, double *buffer)
{
    for(int i = 0; i < length; ++i)
        buffer[i] = randDouble(min[i], max[i]);
}

void fillArrayRandomDoubleSameInterval(double min, double max, int length, double *buffer)
{
    for(int i = 0; i < length; ++i)
        buffer[i] = randDouble(min, max);
}
