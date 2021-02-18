#include <parpecommon/misc.h>
#include <parpecommon/logging.h>

#include <alloca.h>
#include <cassert>
#include <cerrno>
#include <execinfo.h>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <cstring>
#include <sys/stat.h>
#include <sys/types.h>
#include <ctime>
#include <unistd.h>
#include <dlfcn.h> // dladdr
#include <cxxabi.h> // __cxa_demangle
#include <sstream>
#include <cstdlib> // getenv
#include <algorithm>
#include <random>

#include <pthread.h>

#include <gsl/gsl-lite.hpp>

#ifdef PARPE_ENABLE_MPI
#include <mpi.h>
#endif

// void printMatlabArray(const double *buffer, int len)
//{
//    printf("[");
//    printfArray(buffer, len - 1, "%e, ");
//    printf("%e]\n", buffer[len - 1]);
//    fflush(stdout);
//}


namespace parpe {

bool fileExists(const char *name) {
    struct stat buffer{};
    return (stat(name, &buffer) == 0);
}

void createDirectoryIfNotExists(char *dirName) {
    struct stat st = {};

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
    Expects(file_path && *file_path);

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
    Expects(file_path);
    char tmp[strlen(file_path) + 1];

    strncpy(tmp, file_path, sizeof(tmp) -1);
    tmp[sizeof(tmp) - 1] = '\0';

    return mkpath(tmp, mode);
}

void strFormatCurrentLocaltime(gsl::span<char> buffer, const char *format) {
    time_t current_time;
    struct tm local_time;
    time(&current_time);
    localtime_r(&current_time, &local_time);

    strftime(buffer.data(), buffer.size(), format, &local_time);
}

void runInParallelAndWaitForFinish(void *(*function)(void *), void **args,
                                   int numArgs) {
    // create threads
    pthread_attr_t threadAttr;
    pthread_attr_init(&threadAttr);
    pthread_attr_setdetachstate(&threadAttr, PTHREAD_CREATE_JOINABLE);

    auto threads = static_cast<pthread_t *>(alloca(numArgs * sizeof(pthread_t)));

    for (int i = 0; i < numArgs; ++i) {
        pthread_create(&threads[i], &threadAttr, function, args[i]);
    }
    pthread_attr_destroy(&threadAttr);

    // wait for finish
    for (int i = 0; i < numArgs; ++i) {
        pthread_join(threads[i], nullptr);
        logmessage(LOGLVL_DEBUG, "Thread i %d finished", i);
    }
    logmessage(LOGLVL_DEBUG, "All k threads finished.");
}

void printBacktrace(int nMaxFrames) {
    void *array[nMaxFrames];
    size_t size;
    size = backtrace(array, nMaxFrames);
    backtrace_symbols_fd(array, size, STDERR_FILENO);
}

std::string getBacktrace(int nMaxFrames)
{
    std::ostringstream oss;

    void *callstack[nMaxFrames];
    int nFrames = backtrace(callstack, nMaxFrames);
    char **symbols = backtrace_symbols(callstack, nFrames);

    char buf[1024];
    for (int i = 0; i < nFrames; i++) {
        Dl_info info;
        if (dladdr(callstack[i], &info) && info.dli_sname) {
            char *demangled = nullptr;
            int status = -1;
            if (info.dli_sname[0] == '_')
                demangled = abi::__cxa_demangle(info.dli_sname, nullptr, nullptr, &status);
            snprintf(buf, sizeof(buf), "%-3d %*p %s + %zd\n",
                     i, int(2 + sizeof(void*) * 2), callstack[i],
                     status == 0 ? demangled :
                     info.dli_sname == nullptr ? symbols[i] : info.dli_sname,
                     (char *)callstack[i] - (char *)info.dli_saddr);
            free(demangled);
        } else {
            snprintf(buf, sizeof(buf), "%-3d %*p %s\n",
                     i, int(2 + sizeof(void*) * 2), callstack[i], symbols[i]);
        }
        oss << buf;
    }
    free(symbols);
    if (nFrames == nMaxFrames)
        oss << "[truncated]\n";

    return oss.str();

}

double randDouble(double min, double max) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(min, max);
    return dis(gen);
}

void fillArrayRandomDoubleIndividualInterval(gsl::span<const double> min,
                                             gsl::span<const double> max,
                                             gsl::span<double> buffer) {
    Expects(min.size() == max.size());
    Expects(min.size() == buffer.size());

    for (gsl::span<double>::index_type i = 0; i < buffer.size(); ++i)
        buffer[i] = randDouble(min[i], max[i]);
}

void fillArrayRandomDoubleSameInterval(double min, double max,
                                       gsl::span<double> buffer) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(min, max);

    for (gsl::span<double>::index_type i = 0; i < buffer.size(); ++i)
        buffer[i] = dis(gen);
}


int getMpiRank() {
    int mpiRank = -1;
#ifdef PARPE_ENABLE_MPI
    if (getMpiActive()) {
        MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    }
#endif
    return mpiRank;
}

int getMpiCommSize() {
    int mpiCommSize = -1;

#ifdef PARPE_ENABLE_MPI
    if (getMpiActive()) {
        MPI_Comm_size(MPI_COMM_WORLD, &mpiCommSize);
    }
#endif

    return mpiCommSize;
}

int getMpiActive()
{
#ifdef PARPE_ENABLE_MPI
    int result = 0;

    MPI_Initialized(&result);
    if(!result)
        return false;

    MPI_Finalized(&result);
    return !result;
#else
    return false;
#endif
}

void CpuTimer::reset()
{
    start = roundStart = clock();
}

double CpuTimer::getRound()
{
    auto now = clock();
    auto timeRound = static_cast<double>(now - roundStart) / CLOCKS_PER_SEC;
    roundStart = now;

    return timeRound;
}

double CpuTimer::getTotal() const
{
    auto now = clock();
    return (double)(now - start) / CLOCKS_PER_SEC;
}

WallTimer::WallTimer()
{
    reset();
}

void WallTimer::reset()
{
    roundStart = start = std::chrono::system_clock::now();
}

double WallTimer::getRound()
{
    std::chrono::duration<double> duration = (std::chrono::system_clock::now() - roundStart);
    roundStart = std::chrono::system_clock::now();
    return duration.count();
}

double WallTimer::getTotal() const
{
    std::chrono::duration<double> duration = (std::chrono::system_clock::now() - start);
    return duration.count();
}

void finalizeMpiIfNeeded()
{
#ifdef PARPE_ENABLE_MPI
    if(parpe::getMpiActive())
        MPI_Finalize();
#endif
}

bool almostEqual(double a, double b)
{
    if (std::isnan(a) && std::isnan(b))
        return true;

    return std::fabs(a - b) < (std::fabs(a) + std::fabs(b))
            * std::numeric_limits<double>::epsilon();
}


} // namespace parpe
