//
// Created by Mingyu on 23-9-30.
//

#include "utils.h"


bool double_equal(double a, double b, double epsilon) {
    return std::abs(a - b) < epsilon;
}

double get_current_time() {
    timeval t;
    gettimeofday(&t, 0);
    return (double) t.tv_sec + (double) t.tv_usec / 1000000;
}

bool isFileExists_ifstream(char *name) {
    std::ifstream f(name);
    return f.good();
}

bool check_be_cover(unsigned L1, unsigned R1, unsigned L2, unsigned R2) {
    if (L2 <= L1 && R1 <= R2) return true;
    return false;
}
