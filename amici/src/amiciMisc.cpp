#include <logging.h>
#include <sstream>

namespace parpe {

void printAmiciErrMsgIdAndTxt(const char *identifier, const char *format, ...) {
    std::stringstream ss;
    if(identifier) {
        ss <<"["<<identifier<<"] "<<format;
    }

    va_list argptr;
    va_start(argptr,format);
    logmessage(LOGLVL_ERROR, ss.str().c_str(), argptr);
    va_end(argptr);
}

void printAmiciWarnMsgIdAndTxt(const char *identifier, const char *format, ...) {
    std::stringstream ss;
    if(identifier) {
        ss <<"["<<identifier<<"] "<<format;
    }

    va_list argptr;
    va_start(argptr,format);
    logmessage(LOGLVL_WARNING, ss.str().c_str(), argptr);
    va_end(argptr);
}

} // namespace parpe
