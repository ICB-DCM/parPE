#include <parpecommon/logging.h>

#include <parpecommon/parpeConfig.h>
#include <parpecommon/misc.h> // getMpiActive

#include <fstream>
#include <ctime>
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <unistd.h>
#include <sstream>
#include <thread>

#ifdef PARPE_ENABLE_MPI
#include <mpi.h>
#endif

namespace parpe {

const char *loglevelShortStr[] = {"", "CRI", "ERR", "WRN", "INF", "DBG"};
loglevel minimumLogLevel = loglevel::debug;
static void printlogmessage(loglevel lvl, const char *message);

std::string printfToString(const char *fmt, va_list ap) {
    // Get size of string
    va_list ap_count;
    va_copy(ap_count, ap);
    auto size = vsnprintf(nullptr, 0, fmt, ap_count);
    va_end(ap_count);
    ++size;

    // actual formatting
    auto buf = std::make_unique<char []>(size);
    size = vsnprintf(buf.get(), size, fmt, ap);

    return std::string(buf.get(), size);
}

void logmessage(loglevel lvl, const char *format, ...)
{
    va_list argptr;
    va_start(argptr, format);
    auto str = printfToString(format, argptr);
    va_end(argptr);
    logmessage(lvl, str);
}

void logmessage(loglevel lvl, const char *format, va_list argptr) {
    logmessage(lvl, printfToString(format, argptr));
}

void logProcessStats()
{
    std::ifstream file("/proc/self/status");
    if (file.is_open()) {
        std::string line;
        while (std::getline(file, line)) {
            if(line.rfind("Vm", 0) == 0
                    || line.rfind("Rss", 0) == 0) {
                logmessage(loglevel::debug, line);
            }
        }
        file.close();
    }
}

void printMPIInfo() {
#ifdef PARPE_ENABLE_MPI
    int mpiActive = getMpiActive();

    if(mpiActive) {
        int mpiCommSize, mpiRank;
        MPI_Comm_size(MPI_COMM_WORLD, &mpiCommSize);
        MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

        char procName[MPI_MAX_PROCESSOR_NAME];
        int procNameLen;
        MPI_Get_processor_name(procName, &procNameLen);

        logmessage(loglevel::debug, "Rank %d/%d running on %s.", mpiRank,
                   mpiCommSize, procName);
    } else {
        logmessage(loglevel::debug, "MPI not initialized.");
    }
#else
    logmessage(loglevel::debug, "MPI support disabled.");
#endif
}


void printDebugInfoAndWait(int seconds) {
    //int i = 0;
    char hostname[256];
    gethostname(hostname, sizeof(hostname));
    logmessage(loglevel::debug,
               "PID %d on %s ready for attach (will wait for %ds)", getpid(),
               hostname, seconds);
    fflush(stdout);
    //while (0 == i)
    sleep(seconds);
}

void logmessage(loglevel lvl, const std::string &msg)
{
    std::stringstream ss(msg);
    std::string line;

    while(std::getline(ss, line, '\n'))
        printlogmessage(lvl, line.c_str());
}

void printlogmessage(loglevel lvl, const char *message)
{
    if(minimumLogLevel < lvl)
        return;


    // TODO: fileLogLevel, consoleLogLevel
    // Coloring
    switch (lvl) {
    case loglevel::critical:
        printf(ANSI_COLOR_MAGENTA);
        break;
    case loglevel::error:
        printf(ANSI_COLOR_RED);
        break;
    case loglevel::warning:
        printf(ANSI_COLOR_YELLOW);
        break;
    case loglevel::debug:
        printf(ANSI_COLOR_CYAN);
        break;
    case loglevel::info:
        printf(ANSI_COLOR_GREEN);
        break;
    }

    // Timestamp
    char dateBuffer[50];

    strFormatCurrentLocaltime(dateBuffer, "[%Y-%m-%d %H:%M:%S] ");
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
    std::ostringstream thread_id_oss;
    thread_id_oss << std::this_thread::get_id();
    auto thread_id {thread_id_oss.str()};

    printf("[%*d:%s/%s] ", 1 + static_cast<int>(log10(mpiCommSize)),
           mpiRank, thread_id.c_str(), procName );
    printf("%s", message);
    printf("%s\n", ANSI_COLOR_RESET);

    switch (lvl) {
    case loglevel::critical:
        [[fallthrough]];
    case loglevel::error:
        fflush(stdout);
        break;
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
    va_list argptr;
    va_start(argptr, format);
    logmessage(lvl, format, argptr);
    va_end(argptr);
}

void Logger::logmessage(loglevel lvl, const char *format, va_list argptr) const {
    logmessage(lvl, printfToString(format, argptr));
}

void Logger::setPrefix(const std::string &pre) {
    prefix = pre;
}

const std::string &Logger::getPrefix() const {
    return prefix;
}

} // namespace parpe
