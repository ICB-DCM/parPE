#include <parpecommon/logging.h>

#include <parpecommon/parpeConfig.h>
#include <parpecommon/misc.h> // getMpiActive

#include <ctime>
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <unistd.h>
#include <sstream>

#ifdef PARPE_ENABLE_MPI
#include <mpi.h>
#endif

namespace parpe {

const char *loglevelShortStr[] = {"", "CRI", "ERR", "WRN", "INF", "DBG"};
loglevel minimumLogLevel = LOGLVL_DEBUG;

void logmessage(loglevel lvl, const char *format, ...)
{
    va_list argptr;
    va_start(argptr,format);
    logmessage(lvl, format, argptr);
    va_end(argptr);
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
#ifdef PARPE_ENABLE_MPI
    int mpiInitialized = 0;
    MPI_Initialized(&mpiInitialized);

    if(mpiInitialized) {
        int mpiCommSize, mpiRank;
        MPI_Comm_size(MPI_COMM_WORLD, &mpiCommSize);
        MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

        char procName[MPI_MAX_PROCESSOR_NAME];
        int procNameLen;
        MPI_Get_processor_name(procName, &procNameLen);

        logmessage(LOGLVL_DEBUG, "Rank %d/%d running on %s.", mpiRank,
                   mpiCommSize, procName);
    } else {
        logmessage(LOGLVL_DEBUG, "MPI not initialized.");
    }
#else
    logmessage(LOGLVL_DEBUG, "MPI support disabled.");
#endif
}


void printDebugInfoAndWait(int seconds) {
    //int i = 0;
    char hostname[256];
    gethostname(hostname, sizeof(hostname));
    logmessage(LOGLVL_DEBUG,
               "PID %d on %s ready for attach (will wait for %ds)", getpid(),
               hostname, seconds);
    fflush(stdout);
    //while (0 == i)
    sleep(seconds);
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

void logmessage(loglevel lvl, const char *format, va_list argptr)
{
    if(minimumLogLevel < lvl)
        return;


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

    printf("[%s] ", loglevelShortStr[static_cast<int>(lvl)]);

    // MPI info
    int mpiCommSize = 1, mpiRank = -1;
#ifdef PARPE_ENABLE_MPI
    auto mpiActive = getMpiActive();

    if(mpiActive) {
        MPI_Comm_size(MPI_COMM_WORLD, &mpiCommSize);
        MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    }

    char procName[MPI_MAX_PROCESSOR_NAME];
    procName[0] = '\0';
    if(mpiActive) {
        int procNameLen;
        MPI_Get_processor_name(procName, &procNameLen);
    }
#else
    auto procName = "";
#endif
    printf("[%*d/%s] ", 1 + static_cast<int>(log10(mpiCommSize)), mpiRank, procName);

    // Message
    vprintf(format, argptr);
    printf(ANSI_COLOR_RESET "\n");

    switch (lvl) {
    case LOGLVL_CRITICAL:
    case LOGLVL_ERROR:
        fflush(stdout);
    default:
        break;
    }

}

Logger::Logger(std::string prefix) : prefix(std::move(prefix)) {}

std::unique_ptr<Logger> Logger::getChild(const std::string &appendedPrefix) const {
    return std::make_unique<Logger>(prefix + appendedPrefix);
}

void Logger::logmessage(loglevel lvl, const std::string &msg) const {
    parpe::logmessage(lvl, "[" + prefix + "] " + msg);
}

void Logger::logmessage(loglevel lvl, const char *format, ...) const {
    std::string str = "[" + prefix + "] " + format;

    va_list argptr;
    va_start(argptr,format);
    parpe::logmessage(lvl, str.c_str(), argptr);
    va_end(argptr);
}

void Logger::logmessage(loglevel lvl, const char *format, va_list argptr) const {
    std::string str = "[" + prefix + "] " + format;
    parpe::logmessage(lvl, str.c_str(), argptr);
}

void Logger::setPrefix(const std::string &pre) {
    prefix = pre;
}

const std::string &Logger::getPrefix() const {
    return prefix;
}

} // namespace parpe
