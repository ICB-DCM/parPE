#include "testingMisc.h"

int randInt(int min, int max) {
    return min + rand() / (double) RAND_MAX * (max - min);
}
