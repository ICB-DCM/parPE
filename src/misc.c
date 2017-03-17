#include "misc.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <mpi.h>
#include <time.h>
#include <alloca.h>
#include <unistd.h>
#include <math.h>
#include <sys/stat.h>

const char *loglevelShortStr[] = {"", "CRI", "ERR", "WRN", "INF", "DBG"};

void error(const char *message) { // exit?
    logmessage(LOGLVL_ERROR, message);
}

void warning(const char *message) {
    logmessage(LOGLVL_WARNING, message);
}

void getLatinHyperCubeSamples(int numParameters, int numSamples, double *sample) {
    for(int i = 0; i < numParameters; ++i) {
        double tmpSample[numSamples];

        for(int j = 0; j < numSamples; ++j)
            tmpSample[j] = rand();

        double tmpRank[numSamples];

        rank(tmpSample, tmpRank, numSamples);

        for(int j = 0; j < numSamples; ++j) {
            //double add = 0.5; // square center
            double add = rand() / (double)RAND_MAX;
            sample[numParameters * j + i] = (tmpRank[j] + add) / numSamples;
        }
    }
}


int doubleSort(const void *x, const void *y) {
    return (*(double*)x - *(double*)y);
}

void rank(const double *in, double *out, int length) {
    memcpy(out, in, sizeof(double) * length);
    qsort(out, length, sizeof(double), doubleSort);

    for(int i = 0; i < length; ++i) {
        double curVal = out[i];

        for(int j = 0; j < length; ++j) {
            if(in[j] == curVal) {
                out[i] = j;
                break;
            }
        }
    }
}

void logmessage(loglevel lvl, const char *format, ...)
{
    int mpiInitialized = 0;
    MPI_Initialized(&mpiInitialized);

    // TODO: fileLogLevel, consoleLogLevel
    // Coloring
    switch (lvl) {
    case LOGLVL_CRITICAL:
        printf(ANSI_COLOR_MAGENTA);
        break;
    case LOGLVL_ERROR:
        printf(ANSI_COLOR_RED);
        break;
    case LOGLVL_WARNING:
        printf(ANSI_COLOR_YELLOW);
        break;
    case LOGLVL_DEBUG:
        printf(ANSI_COLOR_CYAN);
        break;
    case LOGLVL_INFO:
        printf(ANSI_COLOR_GREEN);
        break;
    }

    // Timestamp
    time_t timer;
    time(&timer);
    struct tm* tm_info;
    tm_info = localtime(&timer);
    char dateBuffer[50];
    strftime(dateBuffer, 25, "[%Y-%m-%d %H:%M:%S] ", tm_info);
    fputs(dateBuffer, stdout);

    printf("[%s] ", loglevelShortStr[lvl]);

    // MPI info
    int mpiCommSize = 1, mpiRank = -1;
    if(mpiInitialized) {
        MPI_Comm_size(MPI_COMM_WORLD, &mpiCommSize);
        MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    }

    char procName[MPI_MAX_PROCESSOR_NAME];
    procName[0] = '\0';
    if(mpiInitialized) {
        int procNameLen;
        MPI_Get_processor_name(procName, &procNameLen);
    }
    printf("[%*d/%s] ", 1 + (int)log10(mpiCommSize), mpiRank, procName);

    // Message
    va_list argptr;
    va_start(argptr,format);
    vprintf(format, argptr);
    va_end(argptr);
    printf(ANSI_COLOR_RESET "\n");

    switch (lvl) {
    case LOGLVL_CRITICAL:
    case LOGLVL_ERROR:
        fflush(stdout);
    default:
        break;
    }
}

//void printMatlabArray(const double *buffer, int len)
//{
//    printf("[");
//    printfArray(buffer, len - 1, "%e, ");
//    printf("%e]\n", buffer[len - 1]);
//    fflush(stdout);
//}

void logProcessStats()
{
    int bufSize = 1024;
    char *buffer = malloc(bufSize);

    FILE* status = fopen( "/proc/self/status", "r" );

    while (fgets(buffer, bufSize, status)) {
        buffer[strlen(buffer) - 1] = '\0'; // remove \n
        logmessage(LOGLVL_DEBUG, "%s", buffer);
    }

    fclose(status);
    free(buffer);
}


int checkGradient(objectiveFunction objFun, objectiveFunctionGradient objFunGrad, int nParams, double *theta, double epsilon, int *indices, int nIndices) {
    double gradient[nParams];
    objFunGrad(theta, gradient);

    double thetaTmp[nParams];
    memcpy(thetaTmp, theta, sizeof(double) * nParams);

    double f = 0;
    objFun(theta, &f);

    printf("Index\tGradient\tfd_f\t\t(delta)\t\tfd_c\t\t(delta)\t\tfd_b\t\t(delta)\n");

    for(int i = 0; i < nIndices; ++i) {
        int curInd = indices[i];
        double fb = 0, ff = 0;

        thetaTmp[curInd] = theta[curInd] + epsilon;
        objFun(thetaTmp, &ff);

        thetaTmp[curInd] = theta[curInd] - epsilon;
        objFun(thetaTmp, &fb);

        double fd_f = (ff - f) / epsilon;

        double fd_b = (f - fb) / epsilon;

        double fd_c = (ff - fb) / (2 * epsilon);

        thetaTmp[curInd] = theta[curInd];

        printf("%d\t%f\t%f\t(%f)\t%f\t(%f)\t%f\t(%f)\n", curInd, gradient[curInd], fd_f, gradient[curInd] - fd_f, fd_c, gradient[curInd] - fd_c, fd_b, gradient[curInd] - fd_b);
    }

    return 0;
}


void printMPIInfo() {
    int mpiCommSize, mpiRank;
    MPI_Comm_size(MPI_COMM_WORLD, &mpiCommSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

    char procName[MPI_MAX_PROCESSOR_NAME];
    int procNameLen;
    MPI_Get_processor_name(procName, &procNameLen);

    logmessage(LOGLVL_DEBUG, "Rank %d/%d running on %s.", mpiRank, mpiCommSize, procName);
}


void printDebugInfoAndWait() {
    //int i = 0;
    char hostname[256];
    gethostname(hostname, sizeof(hostname));
    logmessage(LOGLVL_DEBUG, "PID %d on %s ready for attach", getpid(), hostname);
    fflush(stdout);
    //while (0 == i)
        sleep(15);
}

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

