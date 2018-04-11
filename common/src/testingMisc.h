#ifndef PARPE_TESTING_MISC_H
#define PARPE_TESTING_MISC_H

#include <functional>
#include <string>
#include <iostream>
#include <cstdio>
#include <unistd.h>

namespace parpe {

/**
 * @brief likelihood offset for sigma = 1
 * @param n
 * @return
 */
double getLogLikelihoodOffset(int n);

int randInt(int min, int max);

bool withinTolerance(double expected, double actual, double atol, double rtol, int index);

void checkEqualArray(const double *expected, const double *actual, int length, double atol, double rtol);

std::string captureStreamToString(std::function<void()> f, std::ostream &os = std::cout);

std::string captureStreamToString(std::function<void()> f, std::FILE* captureStream = stdout, int captureStreamFd = STDOUT_FILENO);

} // namespace parpe

#endif
