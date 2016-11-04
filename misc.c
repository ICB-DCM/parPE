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

void getLatinHyperCubeSample(int numParameters, double *sample) {
    double tmpSample[numParameters];

    for(int j = 0; j < numParameters; ++j)
        tmpSample[j] = rand();

    rank(tmpSample, sample, numParameters);

    for(int j = 0; j < numParameters; ++j) {
        // sample[j] = (sample[j] + 0.5) / numParameters; // square center
        sample[j] = (sample[j] + rand() / (double)RAND_MAX) / numParameters;
    }
}

double *getLatinHyperCubeSamples(int numParameters, int numSamples) {
    double *sample = malloc(sizeof(double) * numParameters * numSamples);

    for(int i = 0; i < numSamples; ++i) {
        getLatinHyperCubeSample(numParameters, &sample[i * numParameters]);
    }

    return sample;
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
