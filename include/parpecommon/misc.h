#ifndef CPP_MISC_H
#define CPP_MISC_H

#include <parpecommon/parpeConfig.h>

#include <chrono>
#include <cstdio>
#include <memory>
#include <ostream>
#include <stdlib.h>
#include <vector>

#include <gsl/gsl-lite.hpp>

template<typename T>
std::ostream&
operator<<(std::ostream& o, std::vector<T> const& v)
{
    o << "[ ";
    for (auto const& e : v)
        o << e << " ";
    o << "]";
    return o;
}

template<typename T>
std::ostream&
operator<<(std::ostream& o, gsl::span<T> const& v)
{
    o << "[ ";
    for (auto const& e : v)
        o << e << " ";
    o << "]";
    return o;
}

namespace parpe {

class WallTimer
{
public:
    WallTimer();

    void reset();

    double getRound();

    double getTotal() const;

private:
    std::chrono::time_point<std::chrono::system_clock> start;
    std::chrono::time_point<std::chrono::system_clock> roundStart;
};

class CpuTimer
{
public:
    CpuTimer() = default;

    void reset();

    double getRound();

    double getTotal() const;

private:
    clock_t start = clock();
    clock_t roundStart = clock();
};

#define RELEASE_ASSERT(expr, msg)                                              \
    if (!(expr)) {                                                             \
        /* NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-array-to-pointer-decay, \
         * cppcoreguidelines-pro-type-vararg) */                               \
        printf("CRITICAL: Assertion %s in %s:%d failed (%s)\n",                \
               (#expr),                                                        \
               __FILE__,                                                       \
               __LINE__,                                                       \
               msg);                                                           \
        abort();                                                               \
    }

void
strFormatCurrentLocaltime(gsl::span<char> buffer, const char* format);

void
printBacktrace(int nMaxFrames = 20);

std::string
getBacktrace(int nMaxFrames = 20);

double
randDouble(double min, double max);

/**
 * @brief fillArrayRandomDoubleIndividualInterval Fill "buffer" with
 * random double values, drawn from an interval [min, max] given for each value.
 * @param min
 * @param max
 * @param buffer
 */
void
fillArrayRandomDoubleIndividualInterval(gsl::span<const double> min,
                                        gsl::span<const double> max,
                                        gsl::span<double> buffer);

void
fillArrayRandomDoubleSameInterval(double min,
                                  double max,
                                  gsl::span<double> buffer);

int
getMpiRank();
int
getMpiCommSize();
int
getMpiActive();

void
finalizeMpiIfNeeded();

template<typename T_TEST, typename T_BOUNDS>
bool
withinBounds(long int n,
             T_TEST const* x,
             const T_BOUNDS* min,
             const T_BOUNDS* max)
{
    for (int i = 0; i < n; ++i)
        if (x[i] < min[i])
            return false;

    for (int i = 0; i < n; ++i)
        if (x[i] > max[i])
            return false;

    return true;
}

/**
 * @brief The Like std::unique_lock, but unlocking a mutex on construction and
 * locking on destruction.
 */
template<typename MUTEX>
class InverseUniqueLock
{
  public:
    explicit InverseUniqueLock(MUTEX* mutex)
      : mutex(mutex)
    {
        mutex->unlock();
    }

    InverseUniqueLock(InverseUniqueLock& other) = delete;

    InverseUniqueLock& operator=(const InverseUniqueLock& other) = delete;

    InverseUniqueLock(InverseUniqueLock&& other) noexcept
    {
        mutex = other.mutex;
        other.mutex = nullptr;
    }

    InverseUniqueLock const& operator=(InverseUniqueLock&& fp) = delete;

    ~InverseUniqueLock() { mutex->lock(); }

  private:
    MUTEX* mutex = nullptr;
};

/**
 * @brief Check if a and b are equal to machine precision
 * @param a
 * @param b
 * @return
 */
bool
almostEqual(double a, double b);

} // namespace parpe

#ifndef __cpp_lib_make_unique
// custom make_unique while we are still using c++11
namespace std {
template<typename T, typename... Args>
std::unique_ptr<T>
make_unique(Args&&... args)
{
    return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}
}
#endif

#endif
