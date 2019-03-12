#ifndef CPP_MISC_H
#define CPP_MISC_H

#include <parpecommon/parpeConfig.h>

#include <stdlib.h>
#include <memory>
#include <cstdio>
#include <chrono>
#include <vector>
#include <ostream>

#include <gsl/gsl-lite.hpp>

template<typename T>
std::ostream& operator <<(std::ostream& o, std::vector<T> const& v) {
    o << "[ ";
    for(auto const& e: v)
        o << e << " ";
    o << "]";
    return o;
}

template<typename T>
std::ostream& operator <<(std::ostream& o, gsl::span<T> const& v) {
    o << "[ ";
    for(auto const& e: v)
        o << e << " ";
    o << "]";
    return o;
}

namespace parpe {

class WallTimer {
public:
    WallTimer();

    void reset();

    double getRound();

    double getTotal();

    std::chrono::time_point<std::chrono::system_clock> start;
    std::chrono::time_point<std::chrono::system_clock> roundStart;

};

class CpuTimer {
public:
    CpuTimer() = default;

    void reset();

    double getRound();

    double getTotal();

    clock_t start = clock();
    clock_t roundStart = clock();
};


#define RELEASE_ASSERT(expr, msg) \
    if(!(expr)) { \
        /* NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-array-to-pointer-decay, cppcoreguidelines-pro-type-vararg) */ \
        printf("CRITICAL: Assertion %s in %s:%d failed (%s)\n", \
                          (#expr), __FILE__, __LINE__, msg); \
        abort(); \
    }

// void printMatlabArray(const double *buffer, int len);

/**
 * @brief Check if file or directory exists
 * @param name
 * @return True if exists, false if not
 */
bool fileExists(const char *name);

int mkpath(char *file_path, mode_t mode);

int mkpathConstChar(const char *file_path, mode_t mode);

void createDirectoryIfNotExists(char *dirName);

void strFormatCurrentLocaltime(char *buffer, size_t bufferSize,
                                       const char *format);

void runInParallelAndWaitForFinish(void *(*function)(void *),
                                           void **args, int numArgs);

void printBacktrace(int nMaxFrames = 20);

std::string getBacktrace(int nMaxFrames = 20);

double randDouble(double min, double max);

/**
 * @brief fillArrayRandomDoubleIndividualInterval Fill "buffer" with "length"
 * random double values, drawn from an interval [min, max] given for each value.
 * @param min
 * @param max
 * @param length
 * @param buffer
 */
void fillArrayRandomDoubleIndividualInterval(const double *min,
                                                     const double *max,
                                                     int length,
                                                     double *buffer);

void fillArrayRandomDoubleSameInterval(double min, double max,
                                               int length, double *buffer);

int getMpiRank();
int getMpiCommSize();
int getMpiActive();
bool launchedWithMpi();

void initMpiIfNeeded(int *argc, char ***argv);
void finalizeMpiIfNeeded();

template <typename T_TEST, typename T_BOUNDS>
bool withinBounds(long int n, T_TEST const *x, const T_BOUNDS *min, const T_BOUNDS *max ) {
    for(int i = 0; i < n; ++i)
        if(x[i] < min[i])
            return false;

    for(int i = 0; i < n; ++i)
        if(x[i] > max[i])
            return false;

    return true;
}


/**
 * @brief The Like std::unique_lock, but unlocking a mutex on construction and locking on destruction.
 */
template<typename MUTEX>
class InverseUniqueLock {
public:
    InverseUniqueLock(MUTEX *mutex)
        : mutex(mutex)
    {
        mutex->unlock();
    }

    ~InverseUniqueLock()
    {
        mutex->lock();
    }

private:
    MUTEX *mutex = nullptr;
};

} // namespace parpe

#ifndef __cpp_lib_make_unique
// custom make_unique while we are still using c++11
namespace std {
template<typename T, typename... Args>
std::unique_ptr<T> make_unique(Args&&... args)
{
    return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}
}
#endif

#endif
