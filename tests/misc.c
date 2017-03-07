
double randDouble(double min, double max) {
    return min + rand() / (double) RAND_MAX * (max - min);
}

int randInt(int min, int max) {
    return min + rand() / (double) RAND_MAX * (max - min);
}
