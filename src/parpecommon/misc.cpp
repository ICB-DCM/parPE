#include <parpecommon/misc.h>
#include <parpecommon/logging.h>

#include <execinfo.h>
#include <cmath>
#include <cstdio>
#include <sys/types.h>
#include <ctime>
#include <unistd.h>
#include <dlfcn.h> // dladdr
#include <cxxabi.h> // __cxa_demangle
#include <sstream>
#include <cstdlib> // getenv
#include <algorithm>
#include <random>

#include <gsl/gsl-lite.hpp>

#ifdef PARPE_ENABLE_MPI
#include <mpi.h>
#endif


namespace parpe {

void strFormatCurrentLocaltime(gsl::span<char> buffer, const char *format) {
    time_t current_time;
    struct tm local_time;
    time(&current_time);
    localtime_r(&current_time, &local_time);

    strftime(buffer.data(), buffer.size(), format, &local_time);
}

void printBacktrace(int nMaxFrames) {
    void *array[nMaxFrames];
    auto size = backtrace(array, nMaxFrames);
    backtrace_symbols_fd(array, size, STDERR_FILENO);
}

std::string getBacktrace(int nMaxFrames)
{
    std::ostringstream oss;

    void *callstack[nMaxFrames];
    int nFrames = backtrace(callstack, nMaxFrames);
    auto symbols = std::unique_ptr<char *, void(*)(void*)>{
        backtrace_symbols(callstack, nFrames), free};
    char buf[1024];

    for (int i = 0; i < nFrames; i++) {
        Dl_info info;
        if (dladdr(callstack[i], &info) && info.dli_sname) {
            auto demangled =
                std::unique_ptr<char, void(*)(void*)>{nullptr, free};
            int status = -1;
            if (info.dli_sname[0] == '_')
                demangled.reset(
                    abi::__cxa_demangle(info.dli_sname,nullptr, nullptr,
                                        &status));
            snprintf(buf, sizeof(buf), "%-3d %*p %s + %td\n", i,
                     int(2 + sizeof(void*) * 2), callstack[i],
                     status == 0 ? demangled.get() :
                     info.dli_sname == nullptr ?
                                               symbols.get()[i] : info.dli_sname,
                     (char *)callstack[i] - (char *)info.dli_saddr);
        } else {
            snprintf(buf, sizeof(buf), "%-3d %*p %s\n",
                     i, int(2 + sizeof(void*) * 2), callstack[i], symbols.get()[i]);
        }
        oss << buf;
    }
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

    std::transform(min.begin(), min.end(), max.begin(), buffer.begin(),
                   randDouble);
}

void fillArrayRandomDoubleSameInterval(double min, double max,
                                       gsl::span<double> buffer) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(min, max);

    std::generate(buffer.begin(), buffer.end(), [&]{ return dis(gen); });
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
