#include "logging.h"

#include <time.h>
#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sstream>

namespace parpe {

const char *loglevelShortStr[] = {"", "CRI", "ERR", "WRN", "INF", "DBG"};

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

void logProcessStats()
{
    const int bufSize = 1024;
    char buffer[bufSize];

    FILE* status = fopen( "/proc/self/status", "r" );

    while (fgets(buffer, bufSize, status)) {
        buffer[strlen(buffer) - 1] = '\0'; // remove \n
        logmessage(LOGLVL_DEBUG, "%s", buffer);
    }

    fclose(status);
}

void printMPIInfo() {
    int mpiInitialized = 0;
    MPI_Initialized(&mpiInitialized);

    if(mpiInitialized) {
        int mpiCommSize, mpiRank;
        MPI_Comm_size(MPI_COMM_WORLD, &mpiCommSize);
        MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

        char procName[MPI_MAX_PROCESSOR_NAME];
        int procNameLen;
        MPI_Get_processor_name(procName, &procNameLen);

        logmessage(LOGLVL_DEBUG, "Rank %d/%d running on %s.", mpiRank, mpiCommSize, procName);
    } else {
        logmessage(LOGLVL_DEBUG, "MPI not initialized.");
    }
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

void error(const char *message) { // exit?
    logmessage(LOGLVL_ERROR, message);
}

void warning(const char *message) {
    logmessage(LOGLVL_WARNING, message);
}

void logmessage(loglevel lvl, const std::string &msg)
{
    std::stringstream ss(msg);
    std::string line;

    while(std::getline(ss, line, '\n'))
        logmessage(lvl, line.c_str());
}

} // namespace parpe
