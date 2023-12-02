//
// Created by Mingyu on 23-9-24.
//
#include <chrono>
#include <queue>
#include <unordered_set>
#include <vector>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <unordered_map>
#include <cassert>
#include <chrono>
#include <cstring>
#include <random>
#include <sstream>
#include <boost/dynamic_bitset.hpp>
#include <boost/unordered_map.hpp>
#include <stack>
#include <mutex>
#include <x86intrin.h>
#include <immintrin.h>
#include <malloc.h>
#include <set>
#include <cmath>
#include <queue>
#include <Eigen/Dense>
#include <cassert>
#include <sys/time.h>

#define FLOAT_IVF std::numeric_limits<float>::max()
typedef std::pair<unsigned, unsigned> Segment;
typedef std::lock_guard<std::mutex> LockGuard;
#ifndef SEGMENTSTP_UTILS_H

struct TreeNode {
    std::vector<std::pair<unsigned, double> > nodes_;
    unsigned unique_node;
    int height, fa, depth;
    std::vector<unsigned> pos;
    std::vector<double> dis; //the distance value
    std::vector<unsigned> ch;
};

#define SEGMENTSTP_UTILS_H

bool double_equal(double a, double b, double epsilon = 0.001);

double get_current_time();

bool isFileExists_ifstream(char *name);

bool check_in_segment(unsigned attr, Segment Seg);

bool check_be_cover(unsigned L1,unsigned R1, unsigned L2, unsigned R2);

#endif //SEGMENTSTP_UTILS_H
