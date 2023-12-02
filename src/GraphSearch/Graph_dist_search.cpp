//
// Created by Mingyu on 23-9-25.
//

#include "Graph.h"
#include "utils.h"

using namespace std;

int GraphIndex::bfs_shortest_path(unsigned u, unsigned v, WeightGraph &cur_graph_) {
    std::unordered_map<unsigned, bool> mp;
    std::queue<std::pair<unsigned, int> > Q;
    mp[u] = true;
    Q.emplace(u, 0);
    while (!Q.empty()) {
        unsigned cur_node = Q.front().first;
        int step = Q.front().second;
        if (cur_node == v) return step;
        Q.pop();
        for (auto next_pair: cur_graph_[cur_node]) {
            unsigned next_node = next_pair.first;
            if (mp[next_node]) continue;
            if (next_node == v) {
                return step + 1;
            }
            mp[next_node] = true;
            Q.emplace(next_node, step + 1);
        }
    }
    return -1;
}

double GraphIndex::dijkstra_shortest_path(unsigned u, unsigned v, WeightGraph &cur_graph_) {
    std::vector<bool> mp(cur_graph_.size() + 1, false);
    std::vector<double> dist_mp(cur_graph_.size() + 1, FLOAT_IVF);
    NodeQueue Q;
    Q.emplace(0, u);
    dist_mp[u] = 0;
    while (!Q.empty()) {
        unsigned cur_node = Q.top().second;
        Q.pop();
        if (Q.top().second == v) return dist_mp[v];
        mp[cur_node] = true;
        for (auto next_pair: cur_graph_[cur_node]) {
            unsigned next_node = next_pair.first;
            double next_edge = next_pair.second;
            if (!mp[next_node] && dist_mp[next_node] > dist_mp[cur_node] + next_edge) {
                dist_mp[next_node] = dist_mp[cur_node] + next_edge;
                Q.emplace(dist_mp[next_node], next_node);
            }
        }
    }
    return dist_mp[v];
}

double GraphIndex::bi_dijkstra_shortest_path(unsigned u, unsigned v, WeightGraph &cur_graph_) {
    if (u == v) return 0;
    NodeQueue QF, QB;
    QF.emplace(0, u);
    QB.emplace(0, v);
    std::vector<bool> mpF(cur_graph_.size() + 1, false), mpB(cur_graph_.size() + 1, false);
    std::vector<double> dist_mpF(cur_graph_.size() + 1, FLOAT_IVF);
    std::vector<double> dist_mpB(cur_graph_.size() + 1, FLOAT_IVF);

    dist_mpF[u] = 0;
    dist_mpB[v] = 0;
    unsigned cur_F_node, cur_B_node;
    double cur_F_dist, cur_B_dist;
    double stp_dist = std::numeric_limits<float>::max();
    while (!QF.empty() || !QB.empty()) {
        if (QF.top().first + QB.top().first >= stp_dist) {
            return stp_dist;
        }
        if (!QF.empty()) {
            cur_F_dist = QF.top().first;
            cur_F_node = QF.top().second;
            QF.pop();
            mpF[cur_F_node] = true;
            for (auto next_pair: cur_graph_[cur_F_node]) {
                unsigned next_node = next_pair.first;
                double next_dist = next_pair.second + cur_F_dist;
                if (mpB[next_node] && next_dist + dist_mpB[next_node] < stp_dist) {
                    stp_dist = next_dist + dist_mpB[next_node];
                }
                if (!mpF[next_node]) {
                    if (dist_mpF[next_node] > next_dist) {
                        dist_mpF[next_node] = next_dist;
                        QF.emplace(next_dist, next_node);
                    }
                }
            }
        }
        if (!QB.empty()) {
            cur_B_dist = QB.top().first;
            cur_B_node = QB.top().second;
            QB.pop();
            mpB[cur_B_node] = true;

            for (auto next_pair: cur_graph_[cur_B_node]) {
                unsigned next_node = next_pair.first;
                double next_dist = next_pair.second + cur_B_dist;
                if (mpF[next_node] && next_dist + dist_mpF[next_node] < stp_dist) {
                    stp_dist = next_dist + dist_mpF[next_node];
                }
                if (!mpB[next_node]) {
                    if (dist_mpB[next_node] > next_dist) {
                        dist_mpB[next_node] = next_dist;
                        QB.emplace(next_dist, next_node);
                    }
                }
            }
        }
    }
    return stp_dist;
}

double GraphIndex::PLL_short_dist_query(unsigned u, unsigned v) {
    if (u == v) return 0;
    if (node_order_[u] == -1 || node_order_[v] == -1) return FLOAT_IVF;
    double dist = FLOAT_IVF;
    for (auto cur_pair: Label[u]) {
        unsigned hub = cur_pair.first;
        double hub_dist = cur_pair.second;
        if (Label[v].find(hub) != Label[v].end()) {
            double to_hub_dist = Label[v][hub];
            if (hub_dist + to_hub_dist < dist) {
                dist = hub_dist + to_hub_dist;
            }
        }
    }
    return dist;
}

double GraphIndex::CH_short_dist_query(unsigned u, unsigned v) {
    if (u == v) return 0;
    if (node_order_[u] == -1 || node_order_[v] == -1) return FLOAT_IVF;
    double dist = FLOAT_IVF;
    NodeQueue QF, QB;
    vector<bool> vis_front(node_num_, false);
    vector<bool> vis_back(node_num_, false);
    vector<double> dist_front(node_num_, FLOAT_IVF);
    vector<double> dist_back(node_num_, FLOAT_IVF);
    bool stop_front = false, stop_back = false;
    dist_front[u] = 0;
    QF.emplace(0, u);
    dist_back[v] = 0;
    QB.emplace(0, v);
    while (!QB.empty() || !QF.empty()) {
        if (stop_back && stop_front) break;
        if (stop_front && QB.empty()) break;
        if (stop_back && QF.empty()) break;
        if (!QF.empty() && !stop_front) {
            unsigned cur_node = QF.top().second;
            double cur_dist = QF.top().first;
            QF.pop();
            if (dist_front[cur_node] > dist) stop_front = true;
            vis_front[cur_node] = true;
            if (vis_back[cur_node]) {
                dist = std::min(dist, cur_dist + dist_back[cur_node]);
            }
            for (auto next_pair: CH_graph_[cur_node]) {
                unsigned next_node = next_pair.first;
                double next_dist = next_pair.second + dist_front[cur_node];
                if (!vis_front[next_node]) {
                    if (dist_front[next_node] > next_dist) {
                        dist_front[next_node] = next_dist;
                        QF.emplace(next_dist, next_node);
                    }
                }
            }
        }
        if (!QB.empty() && !stop_back) {
            unsigned cur_node = QB.top().second;
            double cur_dist = QB.top().first;
            QB.pop();
            if (dist_back[cur_node] > dist) stop_back = true;
            vis_back[cur_node] = true;
            if (vis_front[cur_node]) {
                dist = std::min(dist, cur_dist + dist_front[cur_node]);
            }
            for (auto next_pair: CH_graph_[cur_node]) {
                unsigned next_node = next_pair.first;
                double next_dist = next_pair.second + dist_back[cur_node];
                if (!vis_back[next_node]) {
                    if (dist_back[next_node] > next_dist) {
                        dist_back[next_node] = next_dist;
                        QB.emplace(next_dist, next_node);
                    }
                }
            }
        }
    }
    return dist;
}

double GraphIndex::H2H_short_dist_query(unsigned u, unsigned v) {
    if (u == v) return 0;
    if (node_order_[u] == -1 || node_order_[v] == -1) return FLOAT_IVF;
    unsigned rank_u = rank_[u], rank_v = rank_[v];
    unsigned LCA = LCA_query(rank_u, rank_v);
    if (LCA == rank_u)
        return H2H_Tree_[rank_v].dis[H2H_Tree_[rank_u].pos.back()];
    else if (LCA == rank_v)
        return H2H_Tree_[rank_u].dis[H2H_Tree_[rank_v].pos.back()];
    else {
        double tmp = FLOAT_IVF;
        for (int i = 0; i < H2H_Tree_[LCA].pos.size(); i++) {
            if (tmp > H2H_Tree_[rank_u].dis[H2H_Tree_[LCA].pos[i]] + H2H_Tree_[rank_v].dis[H2H_Tree_[LCA].pos[i]])
                tmp = H2H_Tree_[rank_u].dis[H2H_Tree_[LCA].pos[i]] + H2H_Tree_[rank_v].dis[H2H_Tree_[LCA].pos[i]];
        }
        return tmp;
    }
}