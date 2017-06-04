#ifndef MISC_H
#define MISC_H

#ifdef __cplusplus
#define EXTERNC extern "C"
#else
#define EXTERNC
#endif

EXTERNC int randInt(int min, int max);

#endif
