#ifndef MISC_H
#define MISC_H

int randInt(int min, int max);
bool withinTolerance(double expected, double actual, double atol, double rtol, int index);
void checkEqualArray(const double *expected, const double *actual, int length, double atol, double rtol);
#endif
