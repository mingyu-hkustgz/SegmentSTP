//
// Created by mingyu on 23-10-9.
//
#include "SegmentGraph.h"
#include "utils.h"

using namespace std;


/***
 * Note： Only consider undirected graphs
 ***/

void SegmentGraph::add_random_range_attr(unsigned left_range, unsigned right_range) {
    attr_graph_.resize(node_num_);
    srand(123);
    Left_Bound_ = left_range, Right_Bound_ = right_range;
    std::vector<std::unordered_map<unsigned, unsigned> > mp;
    mp.resize(node_num_);
    for (int i = 0; i < node_num_; i++) {
        for (int j = 0; j < weight_graph_[i].size(); j++) {
            unsigned v = weight_graph_[i][j].first;
            if (!mp[i][v]) {
                mp[i][v] = std::rand() % right_range + left_range;
                mp[v][i] = mp[i][v];
            }
            attr_graph_[i].push_back(mp[i][v]);
        }
    }
}
