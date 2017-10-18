#ifndef CPP_MISC_H
#define CPP_MISC_H

#include <stdlib.h>
#include <memory>

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

void shuffle(int *array, size_t numElements);

void runInParallelAndWaitForFinish(void *(*function)(void *),
                                           void **args, int numArgs);

void printBacktrace(int depth);

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

#if __cplusplus < 201402L
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
