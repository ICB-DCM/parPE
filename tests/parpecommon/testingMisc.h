#ifndef PARPE_TESTING_MISC_H
#define PARPE_TESTING_MISC_H

#include <cstdio>
#include <functional>
#include <iostream>
#include <string>
#include <unistd.h>

namespace parpe {

/**
 * @brief likelihood offset for sigma = 1
 * @param n
 * @return
 */
double getLogLikelihoodOffset(int n);

int randInt(int min, int max);

bool withinTolerance(
    double expected,
    double actual,
    double atol,
    double rtol,
    int index);

void checkEqualArray(
    double const* expected,
    double const* actual,
    int length,
    double atol,
    double rtol);

std::string captureStreamToString(
    std::function<void()> const& f,
    std::ostream& os = std::cout);

std::string captureStreamToString(
    std::function<void()> const& f,
    std::FILE* captureStream = stdout,
    int captureStreamFd = STDOUT_FILENO);

} // namespace parpe

#endif
