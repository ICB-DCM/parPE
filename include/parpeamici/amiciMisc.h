#include <amici/defines.h>
#include <amici/misc.h>

namespace parpe {


void printAmiciErrMsgIdAndTxt(const char *identifier, const char *format, ...);

void printAmiciWarnMsgIdAndTxt(const char *identifier, const char *format, ...);

using amici::getUnscaledParameter;

using amici::getScaledParameter;
}
