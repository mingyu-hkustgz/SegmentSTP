//
// Created by Mingyu on 23-10-15.
//
#include "SegmentGraph.h"
#include "Graph.h"
#include "utils.h"
#include "Segment.h"


using namespace std;

bool check_in_segment(unsigned attr, Segment Seg) {
    if (attr < Seg.first) return false;
    if (attr > Seg.second) return false;
    return true;
}

/*
 * use small segment cover big segment
 */
bool Sedge::check_be_cover(Sedge &a) const {
    if (this->L <=a.L && a.R <= this->R) return true;
    return false;
}

bool Sedge::check_in_segment(unsigned left_range, unsigned right_range) const {
    if (left_range <= L && R <= right_range) return true;
    return false;
}



