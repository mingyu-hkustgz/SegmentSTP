//
// Created by Mingyu on 23-9-30.
//
#include "utils.h"
#include "Graph.h"
using namespace std;

void GraphIndex::dijkstra_prune(unsigned u, std::vector<std::pair<unsigned, double>> &vp) {
    NodeQueue Q;
    vector<bool> mp(node_num_, false);
    vector<double> dist_mp(node_num_, FLOAT_IVF);
    dist_mp[u] = 0;
    Q.emplace(0, u);
    while (!Q.empty()) {
        unsigned cur_node = Q.top().second;
        double cur_dist = Q.top().first;
        Q.pop();
        mp[cur_node] = true;
        double temp_dist;
        vector<unsigned> SupNode;
        PLL_SupNode_query(u, cur_node, SupNode, temp_dist);
        if (temp_dist <= cur_dist) continue;
        vp.emplace_back(cur_node, cur_dist);
        for (auto next_pair: weight_graph_[cur_node]) {
            unsigned next_node = next_pair.first;
            double next_dist = next_pair.second + cur_dist;
            if (!mp[next_node]) {
                if (dist_mp[next_node] > next_dist) {
                    dist_mp[next_node] = next_dist;
                    Q.emplace(next_dist, next_node);
                }
            }
        }
    }
}


void GraphIndex::PLL_construct() {
    sort_node_by_degree();
    Label.resize(node_num_);
    unsigned id;
    int cnt = 0;
    /*
     * use degree decent order
     */
    for (int i = (int) node_num_-1; i >= 0; i--) {
        id = invert_node_order_[i];
        if (i % 100 == 0) std::cerr << "PLL current:: " << i << endl;
        vector<pair<unsigned, double>> vp;
        dijkstra_prune(id, vp);
        cnt += 1;
        for (auto &j: vp) {
            Label[j.first].emplace(id, j.second);
        }
    }
}

double GraphIndex::PLL_SupNode_query(unsigned u, unsigned v, vector<unsigned> &SupNode, double &dist) {
    dist = FLOAT_IVF;
    for (auto cur_pair: Label[u]) {
        unsigned hub = cur_pair.first;
        double hub_dist = cur_pair.second;
        if (Label[v].find(hub) != Label[v].end()) {
            double to_hub_dist = Label[v][hub];
            if (hub_dist + to_hub_dist < dist) {
                dist = hub_dist + to_hub_dist;
                SupNode.clear();
                SupNode.push_back(hub);
            } else if (double_equal(dist, hub_dist + to_hub_dist)) {
                SupNode.push_back(hub);
            }
        }
    }
    return dist;
}

