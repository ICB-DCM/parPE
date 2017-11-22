#include "misc.h"
#include "logging.h"
#include <alloca.h>
#include <assert.h>
#include <errno.h>
#include <execinfo.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>
#include <mpi.h>

// void printMatlabArray(const double *buffer, int len)
//{
//    printf("[");
//    printfArray(buffer, len - 1, "%e, ");
//    printf("%e]\n", buffer[len - 1]);
//    fflush(stdout);
//}

namespace parpe {


bool fileExists(const char *name) {
    struct stat buffer;
    return (stat(name, &buffer) == 0);
}

void createDirectoryIfNotExists(char *dirName) {
    struct stat st = {0};

    if (stat(dirName, &st) == -1) {
        mkdir(dirName, 0700);
    }
}

/**
 * @brief Create the path to the given file. Does not error if path exists.
 * @param file_path File name
 * @param mode File mode
 * @return 0 on success, -1 otherwise
 */
int mkpath(char *file_path, mode_t mode) {
    assert(file_path && *file_path);

    for (char *p = strchr(file_path + 1, '/'); p; p = strchr(p + 1, '/')) {
        *p = '\0';
        if (mkdir(file_path, mode) == -1) {
            if (errno != EEXIST) {
                *p = '/';
                return -1;
            }
        }
        *p = '/';
    }
    return 0;
}

int mkpathConstChar(const char *file_path, mode_t mode) {
    assert(file_path && *file_path);
    char tmp[strlen(file_path) + 1];
    strcpy(tmp, file_path);
    return mkpath(tmp, mode);
}

void strFormatCurrentLocaltime(char *buffer, size_t bufferSize,
                               const char *format) {
    time_t timer;
    time(&timer);

    struct tm *tm_info;
    tm_info = localtime(&timer);

    strftime(buffer, bufferSize, format, tm_info);
}

#include <pthread.h>

void runInParallelAndWaitForFinish(void *(*function)(void *), void **args,
                                   int numArgs) {
    // create threads
    pthread_attr_t threadAttr;
    pthread_attr_init(&threadAttr);
    pthread_attr_setdetachstate(&threadAttr, PTHREAD_CREATE_JOINABLE);

    pthread_t *threads = (pthread_t *) alloca(numArgs * sizeof(pthread_t));

    for (int i = 0; i < numArgs; ++i) {
        pthread_create(&threads[i], &threadAttr, function, args[i]);
    }
    pthread_attr_destroy(&threadAttr);

    // wait for finish
    for (int i = 0; i < numArgs; ++i) {
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
    return min + rand() / (double)RAND_MAX * (max - min);
}

void fillArrayRandomDoubleIndividualInterval(const double *min,
                                             const double *max, int length,
                                             double *buffer) {
    for (int i = 0; i < length; ++i)
        buffer[i] = randDouble(min[i], max[i]);
}

void fillArrayRandomDoubleSameInterval(double min, double max, int length,
                                       double *buffer) {
    for (int i = 0; i < length; ++i)
        buffer[i] = randDouble(min, max);
}


int getMpiRank() {
    int mpiRank = -1;

    if (getMpiActive()) {
        MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    }

    return mpiRank;
}

int getMpiCommSize() {
    int mpiCommSize = -1;

    if (getMpiActive()) {
        MPI_Comm_size(MPI_COMM_WORLD, &mpiCommSize);
    }

    return mpiCommSize;
}

int getMpiActive()
{
    int result = 0;

    MPI_Initialized(&result);
    if(!result)
        return false;

    MPI_Finalized(&result);
    return !result;
}

} // namespace parpe
