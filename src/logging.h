#ifndef LOGGING_H
#define LOGGING_H

#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_RESET   "\x1b[0m"

typedef enum loglevel_tag {LOGLVL_CRITICAL = 1, LOGLVL_ERROR, LOGLVL_WARNING, LOGLVL_INFO, LOGLVL_DEBUG} loglevel;

void logmessage(loglevel lvl, const char *format, ...);

void logProcessStats();

void printMPIInfo();

void printDebugInfoAndWait();

// TODO remove
void error(const char *message);
// TODO remove
void warning(const char *message);

#endif
