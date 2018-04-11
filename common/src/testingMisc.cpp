#include "testingMisc.h"
#include <stdlib.h>
#include <cstdio>
#include <cmath>
#include "CppUTest/TestHarness.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <streambuf>
#include <fcntl.h> // O_WRONLY
#include <cassert>

namespace parpe {

bool withinTolerance(double expected, double actual, double atol, double rtol, int index) {
    bool withinTol =  fabs(expected - actual) <= atol || fabs((expected - actual) / (rtol + expected)) <= rtol;

    if(!withinTol && std::isnan(expected) && std::isnan(actual))
        withinTol = true;

    if(!withinTol && std::isinf(expected) && std::isinf(actual))
        withinTol = true;

    if(!withinTol) {
        fprintf(stderr, "ERROR: Expected value %e, but was %e at index %d.\n",expected, actual, index);
        fprintf(stderr, "       Relative error: %e (tolerance was %e)\n", fabs((expected - actual) / (rtol + expected)), rtol);
        fprintf(stderr, "       Absolute error: %e (tolerance was %e)\n", fabs(expected - actual), atol);
        //printBacktrace(12);
    }

    return withinTol;
}

void checkEqualArray(const double *expected, const double *actual, int length, double atol, double rtol) {
    if(!expected && !actual)
        return;

    CHECK_TRUE(expected && actual);

    for(int i = 0; i < length; ++i)
    {
        bool withinTol = withinTolerance(expected[i], actual[i], atol, rtol, i);
        CHECK_TRUE(withinTol);
    }
}

int randInt(int min, int max) {
    return min + rand() / (double) RAND_MAX * (max - min);
}

std::string captureStreamToString(std::function<void()> f, std::ostream &os) {
    std::streambuf* oldOStreamBuf = os.rdbuf();
    os.flush();

    std::ostringstream strOs;
    os.rdbuf(strOs.rdbuf());

    f();

    strOs.flush();
    os.rdbuf( oldOStreamBuf );

    return strOs.str();
}

std::string captureStreamToString(std::function<void()> f, std::FILE* captureStream, int captureStreamFd) {
    // use fmemopen instead of file?

    char tempFileName[TMP_MAX];
    std::tmpnam(tempFileName);

    int newStreamFd = open(tempFileName, O_CREAT | O_WRONLY, S_IRWXU);
    assert(newStreamFd >= 0);

    int oldStreamFd = dup(captureStreamFd);
    assert(oldStreamFd >= 0);
    fflush(captureStream);

    dup2(newStreamFd, captureStreamFd); // replace original fd by tmp file
    close(newStreamFd);  // close remaining copy

    f();
    fflush(captureStream);

    dup2(oldStreamFd, captureStreamFd); // restore (closes tmp file)
    close(oldStreamFd); // close remainingv copy

    std::ifstream ifs(tempFileName);

    return std::string ((std::istreambuf_iterator<char>(ifs)),
                     std::istreambuf_iterator<char>());
}

double getLogLikelihoodOffset(int n) {
    const double pi = atan(1) * 4.0;
    return - n * 0.5 * log(2.0 * pi);
}

} // namespace parpe
