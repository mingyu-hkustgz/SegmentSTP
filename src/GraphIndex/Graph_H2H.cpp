//
// Created by Mingyu on 23-10-1.
//
#include "utils.h"
#include "Graph.h"

using namespace std;


int GraphIndex::match(vector<std::pair<unsigned, double>> &nodes) {
    unsigned nearest = nodes[0].first;
    for (int i = 1; i < nodes.size(); i++) {
        if (rank_[nodes[i].first] > rank_[nearest])
            nearest = nodes[i].first;
    }
    int fa = (int) rank_[nearest];
    return fa;
}

void GraphIndex::build_RMQ_dfs_order(unsigned u) {
    toRMQ[u] = EulerSeq.size();
    EulerSeq.push_back(u);
    for (auto v: H2H_Tree_[u].ch) {
        build_RMQ_dfs_order(v);
        EulerSeq.push_back(u);
    }
}

void GraphIndex::build_RMQ() {
    toRMQ.resize(node_num_);
    build_RMQ_dfs_order(0);
    RMQIndex.push_back(EulerSeq);
    unsigned m = EulerSeq.size();
    for (unsigned i = 2, k = 1; i < m; i = i * 2, k++) {
        vector<unsigned> tmp;
        tmp.resize(m);
        for (int j = 0; j < m - i; j++) {
            unsigned x = RMQIndex[k - 1][j], y = RMQIndex[k - 1][j + i / 2];
            if (H2H_Tree_[x].height < H2H_Tree_[y].height)
                tmp[j] = x;
            else tmp[j] = y;
        }
        RMQIndex.push_back(tmp);
    }
}


unsigned GraphIndex::LCA_query(unsigned u, unsigned v) {
    u = toRMQ[u];
    v = toRMQ[v];
    if (u > v) std::swap(u, v);
    unsigned len = (v - u + 1);
    int i = 1, k = 0;
    while (i * 2 < len) i *= 2, k++;
    v = v - i + 1;
    if (H2H_Tree_[RMQIndex[k][u]].height < H2H_Tree_[RMQIndex[k][v]].height)
        return RMQIndex[k][u];
    else return RMQIndex[k][v];

}


void GraphIndex::build_H2H_index() {
    if (!is_load_ch) {
        CH_decomposition();
        is_load_ch = true;
    }
    std::cerr << "finished CH decomposition" << std::endl;
    build_tree();
    std::cerr << "finished tree decomposition" << std::endl;
    build_tree_index();
    std::cerr << "finished tree index" << std::endl;
}


void GraphIndex::build_tree() {
    node_to_tree_map_.resize(node_num_);
    rank_.resize(node_num_);
    height_max_ = 0;
    int node_len = (int) invert_node_order_.size();
    node_len--;
    while (invert_node_order_[node_len] == -1) {
        node_len--;
    }
    unsigned cur_node = invert_node_order_[node_len];
    TreeNode cur_root{CH_graph_[cur_node], cur_node, 1, -1};
    rank_[cur_node] = 0;
    H2H_Tree_.push_back(cur_root);
    node_len--;
    for (; node_len >= 0; node_len--) {
        cur_node = invert_node_order_[node_len];
        int father = match(CH_graph_[cur_node]);
        H2H_Tree_[father].ch.push_back(H2H_Tree_.size());
        TreeNode node{CH_graph_[cur_node], cur_node, H2H_Tree_[father].height + 1, father, H2H_Tree_[father].depth + 1};
        for (auto next_pair: CH_graph_[cur_node]) {
            unsigned next_node = next_pair.first;
            node_to_tree_map_[next_node].push_back(H2H_Tree_.size());
            if (H2H_Tree_[rank_[next_node]].depth < H2H_Tree_[father].height + 1)
                H2H_Tree_[rank_[next_node]].depth = H2H_Tree_[father].height + 1;
        }
        if (node.height > height_max_) height_max_ = node.height;
        rank_[cur_node] = H2H_Tree_.size();
        H2H_Tree_.push_back(node);
    }
}


void GraphIndex::build_tree_index_dfs(unsigned u, vector<unsigned> &list) {
    unsigned neighbor_size = H2H_Tree_[u].nodes_.size();
    H2H_Tree_[u].pos.resize(neighbor_size + 1);
    H2H_Tree_[u].dis.resize(list.size(), FLOAT_IVF);

    for (unsigned i = 0; i < neighbor_size; i++) {
        for (unsigned j = 0; j < list.size(); j++) {
            if (H2H_Tree_[u].nodes_[i].first == list[j]) {
                H2H_Tree_[u].pos[i] = j;
                H2H_Tree_[u].dis[j] = H2H_Tree_[u].nodes_[i].second;
                break;
            }
        }
    }
    H2H_Tree_[u].pos[neighbor_size] = list.size();


    for (unsigned i = 0; i < neighbor_size; i++) {
        unsigned cur_node = H2H_Tree_[u].nodes_[i].first;
        double cur_dist = H2H_Tree_[u].nodes_[i].second;
        unsigned k = H2H_Tree_[u].pos[i];
        for (unsigned j = 0; j < list.size(); j++) {
            unsigned next_node = list[j];
            double dist_to_next;
            if (k != j) {
                if (k < j)
                    dist_to_next = H2H_Tree_[rank_[next_node]].dis[k];
                else
                    dist_to_next = H2H_Tree_[rank_[cur_node]].dis[j];

                if (H2H_Tree_[u].dis[j] > dist_to_next + cur_dist) {
                    H2H_Tree_[u].dis[j] = dist_to_next + cur_dist;
                }
            }
        }
    }

    list.push_back(H2H_Tree_[u].unique_node);
    for (auto v: H2H_Tree_[u].ch) {
        build_tree_index_dfs(v, list);
    }
    list.pop_back();
}


void GraphIndex::build_tree_index() {
    build_RMQ();
    vector<unsigned> list;
    list.push_back(H2H_Tree_[0].unique_node);
    H2H_Tree_[0].pos.clear();
    H2H_Tree_[0].pos.push_back(0);

    for (auto u: H2H_Tree_[0].ch) {
        build_tree_index_dfs(u, list);
    }

}

