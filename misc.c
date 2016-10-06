#include "misc.h"
#include <stdio.h>

void error(const char *message) { // exit?
    printf("ERROR: %s\n", message);
}

void warning(const char *message) {
    printf("WARNING: %s\n", message);
}
