#ifndef LOGGING_H
#define LOGGING_H

#include <string>
#include <cstdarg>
#include <memory>

namespace parpe {

constexpr const char ANSI_COLOR_RED[] = "\x1b[31m";
constexpr const char ANSI_COLOR_GREEN[] = "\x1b[32m";
constexpr const char ANSI_COLOR_YELLOW[] = "\x1b[33m";
constexpr const char ANSI_COLOR_BLUE[] = "\x1b[34m";
constexpr const char ANSI_COLOR_MAGENTA[] = "\x1b[35m";
constexpr const char ANSI_COLOR_CYAN[] = "\x1b[36m";
constexpr const char ANSI_COLOR_RESET[] = "\x1b[0m";

std::string printfToString(const char *fmt, va_list ap);

enum class loglevel {
    critical = 1,
    error,
    warning,
    info,
    debug
};

// Minimum log level that will be printed
extern loglevel minimumLogLevel;

void logmessage(loglevel lvl, std::string const& msg);
void logmessage(loglevel lvl, const char *format, ...);
void logmessage(loglevel lvl, const char *format, va_list argptr);

/**
 * @brief Print process statistics from /proc/self/status
 */

void logProcessStats();

void printMPIInfo();

void printDebugInfoAndWait(int seconds = 15);

// TODO remove
void error(const char *message);
// TODO remove
void warning(const char *message);

class Logger {
public:
    Logger() = default;
    explicit Logger(std::string prefix);

    std::unique_ptr<Logger> getChild(std::string const& appendedPrefix) const;

    // TODO add stream operator

    void logmessage(loglevel lvl, std::string const& msg) const;
    void logmessage(loglevel lvl, const char *format, ...) const;
    void logmessage(loglevel lvl, const char *format, va_list argptr) const;
    void setPrefix(std::string const& pre);
    std::string const& getPrefix() const;

private:
    std::string prefix;
};

} // namespace parpe
#endif
