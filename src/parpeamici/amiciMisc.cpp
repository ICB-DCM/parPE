#include <parpeamici/amiciMisc.h>

#include <parpecommon/logging.h>

#include <sstream>

namespace parpe {

void printAmiciErrMsgIdAndTxt(const char *identifier, const char *format, ...) {
    std::stringstream ss;
    if(identifier) {
        ss <<"["<<identifier<<"] "<<format;
    }

    va_list argptr;
    va_start(argptr, format); // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay, cppcoreguidelines-pro-type-vararg)
    logmessage(LOGLVL_ERROR, ss.str().c_str(), argptr);
    va_end(argptr); // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay)
}

void printAmiciWarnMsgIdAndTxt(const char *identifier, const char *format, ...) {
    std::stringstream ss;
    if(identifier) {
        ss <<"["<<identifier<<"] "<<format;
    }

    va_list argptr;
    va_start(argptr, format); // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay, cppcoreguidelines-pro-type-vararg)
    logmessage(LOGLVL_WARNING, ss.str().c_str(), argptr);
    va_end(argptr); // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay)
}

} // namespace parpe
