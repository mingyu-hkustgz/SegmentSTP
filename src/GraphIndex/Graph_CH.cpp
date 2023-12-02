//
// Created by Mingyu on 23-9-27.
//
#include "utils.h"
#include "Graph.h"

using namespace std;

void GraphIndex::eliminate_E_order(unsigned u, unsigned v) {
    /*
     * bilateral
     */
    if (E[u].find(v) != E[u].end()) {
        E[u].erase(E[u].find(v));
    }
    if (E[v].find(u) != E[v].end()) {
        E[v].erase(E[v].find(u));
    }
}

void GraphIndex::insert_E_order(unsigned u, unsigned v, double w) {
    if (E[u].find(v) == E[u].end()) {
        E[u].insert(make_pair(v, w));
    } else {
        if (E[u][v] > w)
            E[u][v] = w;
    }
}

void GraphIndex::CH_decomposition() {
    sort_node_by_degree();
    CH_graph_.resize(node_num_);
    CH_short_cut_.resize(node_num_);
    E.resize(node_num_);
    for (int i = 0; i < weight_graph_.size(); i++) {
        for (auto u: weight_graph_[i]) {
            E[i].insert(make_pair(u.first, u.second));
        }
    }
    vector<bool> visited(node_num_, false);
    for (int order = 0; order < node_num_; order++) {
        if (order % 100 == 0) std::cerr << "current order:: " << order << endl;
        int cur_node = invert_node_order_[order];
        if(cur_node == -1) continue;
        visited[cur_node] = true;
        vector<pair<unsigned, double> > cur_neighbor;
        for (auto u: E[cur_node]) {
            if (!visited[u.first]) {
                cur_neighbor.emplace_back(u);
            }
        }
        CH_graph_[cur_node] = cur_neighbor;
        for (auto &u: cur_neighbor) {
            eliminate_E_order(cur_node, u.first);
        }
        if (cur_neighbor.size() < 100) {
            for (int i = 0; i < cur_neighbor.size(); i++) {
                unsigned id_i = cur_neighbor[i].first;
                double weight_i = cur_neighbor[i].second;
                for (auto &j: cur_neighbor) {
                    unsigned id_j = j.first;
                    double weight_j = j.second;
                    insert_E_order(id_i, id_j, weight_i + weight_j);
                    if (id_i < id_j)
                        CH_short_cut_[id_i][id_j].push_back(cur_node);
                }
            }
        } else {
#pragma omp parallel for
            for (int i = 0; i < cur_neighbor.size(); i++) {
                unsigned id_i = cur_neighbor[i].first;
                double weight_i = cur_neighbor[i].second;
                for (auto &j: cur_neighbor) {
                    unsigned id_j = j.first;
                    double weight_j = j.second;
                    insert_E_order(id_i, id_j, weight_i + weight_j);
                    if (id_i < id_j)
                        CH_short_cut_[id_i][id_j].push_back(cur_node);
                }
            }
        }
    }
}