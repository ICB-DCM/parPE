#include "misc.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <include/symbolic_functions.h>

void error(const char *message) { // exit?
    printf("ERROR: %s\n", message);
}

void warning(const char *message) {
    printf("WARNING: %s\n", message);
}

void getLatinHyperCubeSamples(int numParameters, int numSamples, double *sample) {
    for(int i = 0; i < numParameters; ++i) {
        double tmpSample[numSamples];

        for(int j = 0; j < numSamples; ++j)
            tmpSample[j] = rand();

        double tmpRank[numSamples];

        rank(tmpSample, tmpRank, numSamples);

        for(int j = 0; j < numSamples; ++j) {
            //double add = 0.5; // square center
            double add = rand() / (double)RAND_MAX;
            sample[numParameters * j + i] = (tmpRank[j] + add) / numSamples;
        }
    }
}


int doubleSort(const void *x, const void *y) {
    return (*(double*)x - *(double*)y);
}

void rank(const double *in, double *out, int length) {
    memcpy(out, in, sizeof(double) * length);
    qsort(out, length, sizeof(double), doubleSort);

    for(int i = 0; i < length; ++i) {
        double curVal = out[i];

        for(int j = 0; j < length; ++j) {
            if(in[j] == curVal) {
                out[i] = j;
                break;
            }
        }
    }
}
