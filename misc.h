#ifndef CPP_MISC_H
#define CPP_MISC_H

#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_RESET   "\x1b[0m"

void error(const char *message);

void warning(const char *message);

void getLatinHyperCubeSamples(int numParameters, int numSamples, double *sample);

int doubleSort(const void *x, const void *y);

void rank(const double *in, double *out, int length);

typedef enum loglevel_tag {LOGLVL_CRITICAL = 1, LOGLVL_ERROR, LOGLVL_WARNING, LOGLVL_INFO, LOGLVL_DEBUG} loglevel;

void logmessage(loglevel lvl, const char *format, ...);

void printMatlabArray(const double *buffer, int len);
#endif
