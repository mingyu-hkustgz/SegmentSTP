//
// Created by Mingyu on 23-9-27.
//
#include "utils.h"
#include "Graph.h"
#include "SegmentGraph.h"

using namespace std;

void GraphIndex::sort_node_by_degree() {
    std::vector<unsigned> node_degree(node_num_);
    PairQueue Q;
    node_order_.resize(node_num_);
    invert_node_order_.resize(node_num_);
    for (int i = 0; i < node_num_; i++) {
        node_degree[i] = weight_graph_[i].size();
        if (node_degree[i] == 0) node_order_[i] = -1;
        else Q.emplace(node_degree[i], i);
        invert_node_order_[i] = -1;
    }
    int cnt = 0;
    while (!Q.empty()) {
        unsigned id = Q.top().second;
        Q.pop();
        node_order_[id] = cnt;
        invert_node_order_[cnt] = (int) id;
        cnt++;
    }
}

vector<double> SegmentGraph::tree_statistics_analysis() {
    vector<double> res;
    double max_tree_height = 0, sum_tree_height = 0;
    double max_tree_width = 0, sum_tree_width = 0;
    double max_seg_size = 0, sum_seg_size = 0, count_seg_item;
    queue<unsigned> Q;
    for(auto v:roots){
        Q.push(v);
    }

    while(!Q.empty())
    {
        auto v = Q.front();
        Q.pop();

        sum_tree_width += (double) tree_nodes_[v].edges.size() + 1;
        max_tree_width = max(max_tree_width, (double) tree_nodes_[v].edges.size() + 1.0);
        sum_tree_height += tree_nodes_[v].height;
        max_tree_height = max(max_tree_width, (double)tree_nodes_[v].height);
        for(const auto& seg:tree_nodes_[v].edges){
            max_seg_size = max(max_seg_size, (double) seg.second.size());
            count_seg_item += 1.0;
            sum_seg_size += (double) seg.second.size();
        }
        for(auto son:tree_nodes_[v].ch){
            Q.push(son);
        }
    }
    res.push_back(max_tree_width);
    res.push_back(sum_tree_width);
    res.push_back(max_tree_height);
    res.push_back(sum_tree_height);
    res.push_back(max_seg_size);
    res.push_back(sum_seg_size/count_seg_item);
    return res;
}


