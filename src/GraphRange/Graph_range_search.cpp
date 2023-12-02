//
// Created by Mingyu on 23-10-9.
//

#include "utils.h"
#include "SegmentGraph.h"

using namespace std;


double SegmentGraph::dijkstra_range_path(unsigned u, unsigned v, WeightGraph &cur_graph_, Segment seg) {
    std::vector<bool> mp(cur_graph_.size() + 1, false);
    std::vector<double> dist_mp(cur_graph_.size() + 1, FLOAT_IVF);
    NodeQueue Q;
    Q.emplace(0, u);
    dist_mp[u] = 0;
    while (!Q.empty()) {
        unsigned cur_node = Q.top().second;
        Q.pop();
        if (cur_node == v) return dist_mp[v];
        mp[cur_node] = true;
        for (int j = 0; j < cur_graph_[cur_node].size(); j++) {
            if (!check_in_segment(attr_graph_[cur_node][j], seg)) continue;
            auto next_pair = cur_graph_[cur_node][j];
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

double SegmentGraph::bi_dijkstra_range_path(unsigned u, unsigned v, WeightGraph &cur_graph_, Segment seg) {
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
    double stp_dist = FLOAT_IVF;
    while (!QF.empty() || !QB.empty()) {
        if (QF.top().first + QB.top().first >= stp_dist) {
            return stp_dist;
        }
        if (!QF.empty()) {
            cur_F_dist = QF.top().first;
            cur_F_node = QF.top().second;
            QF.pop();
            mpF[cur_F_node] = true;
            for (int j = 0; j < cur_graph_[cur_F_node].size(); j++) {
                if (!check_in_segment(attr_graph_[cur_F_node][j], seg)) continue;
                auto next_pair = cur_graph_[cur_F_node][j];
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
            for (int j = 0; j < cur_graph_[cur_B_node].size(); j++) {
                if (!check_in_segment(attr_graph_[cur_B_node][j], seg)) continue;
                auto next_pair = cur_graph_[cur_B_node][j];
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


double SegmentGraph::range_PLL_query(unsigned u, unsigned v, unsigned left_range, unsigned right_range) {
    if (u == v) return 0;
    double dist = FLOAT_IVF;
    auto s = hopIndex[u];
    auto t = hopIndex[v];
    auto sit = s.begin(), dit = t.begin();
    while (sit != s.end() && dit != t.end()) {
        if (sit->first < dit->first)++sit;
        else if (dit->first < sit->first)++dit;
        else {//same vertex
            for (auto &sp: sit->second.list) {
                for (auto &dp: dit->second.list) {
                    if (check_be_cover(std::min(sp.L, dp.L), std::max(sp.R, dp.R), left_range, right_range)) {
                        dist = std::min(dist, sp.weight + dp.weight);
                    }
                }
            }
            ++sit, ++dit;
        }
    }
    return dist;
}

pair<unsigned int, unsigned int> SegmentGraph::qLCA(unsigned s, unsigned t) {
    unsigned sAncestor = s, tAncestor = t;
    while (true) {
        if (tree_nodes_[sAncestor].height > tree_nodes_[tAncestor].height) {
            if (static_cast<int>(tAncestor) == tree_nodes_[sAncestor].father) return make_pair(tAncestor, tAncestor);
            sAncestor = static_cast<unsigned int>(tree_nodes_[sAncestor].father);
        } else if (tree_nodes_[sAncestor].height < tree_nodes_[tAncestor].height) {
            if (static_cast<int>(sAncestor) == tree_nodes_[tAncestor].father) return make_pair(sAncestor, sAncestor);
            tAncestor = static_cast<unsigned int>(tree_nodes_[tAncestor].father);
        } else {
            if (tree_nodes_[sAncestor].father == -1 || sAncestor == tAncestor)return make_pair(sAncestor, tAncestor);
            sAncestor = static_cast<unsigned int>(tree_nodes_[sAncestor].father);
            tAncestor = static_cast<unsigned int>(tree_nodes_[tAncestor].father);
        }
    }
}


vector<pair<unsigned, double> >
SegmentGraph::down_top_search(unsigned down, unsigned up, unsigned left_range, unsigned right_range) {
    vector<pair<unsigned, double>> result;
    unordered_map<unsigned, double> distmap;

    distmap[down] = 0;
    auto anc = down;
    while (anc != up) {
        unordered_map<unsigned int, double> tmp;
        SegAttr &edges = tree_nodes_[anc].edges;
        double dist = FLOAT_IVF;
        if (distmap.find(anc) != distmap.end())
            dist = distmap[anc];
        for (auto &u: edges) {
            for (auto &medge: u.second.list) {
                if (medge.check_in_segment(left_range, right_range)) {
                    tmp.insert(make_pair(u.first, medge.weight + dist));
                    break;
                }
            }
            if (distmap.find(u.first) != distmap.end()) {
                if (tmp.find(u.first) != distmap.end())
                    tmp[u.first] = min(tmp[u.first], distmap[u.first]);
                else tmp[u.first] = distmap[u.first];
            }
        }
        distmap = std::move(tmp);
        anc = static_cast<unsigned>(tree_nodes_[anc].father);
    }
    for (auto &p: distmap)
        result.emplace_back(p);
    return result;
}


double SegmentGraph::range_tree_decompose_query(unsigned u, unsigned v, unsigned left_range, unsigned right_range) {
    pair<unsigned, unsigned> lca = qLCA(u, v);
    vector<pair<unsigned int, double>> s_result, t_result;
    s_result = down_top_search(u, lca.first, left_range, right_range);
    t_result = down_top_search(v, lca.second, left_range, right_range);
    double distance = FLOAT_IVF;

    for (auto &a: s_result) {
        if (a.second == FLOAT_IVF) continue;
        for (auto &b: t_result) {
            if (b.second == FLOAT_IVF) continue;
            if (a.first == b.first) distance = min(distance, a.second + b.second);
            else {
                distance = min(distance, a.second + b.second + cp_distv1(a.first, b.first, left_range, right_range));
            }
        }
    }
    return distance;
}


pair<unsigned int, unsigned int> SegmentGraph::brute_qLCA(unsigned s, unsigned t, unsigned id) {
    unsigned sAncestor = s, tAncestor = t;
    while (true) {
        if (brute_tree_nodes[id][sAncestor].height > brute_tree_nodes[id][tAncestor].height) {
            if (static_cast<int>(tAncestor) == brute_tree_nodes[id][sAncestor].father)
                return make_pair(tAncestor, tAncestor);
            sAncestor = static_cast<unsigned int>(brute_tree_nodes[id][sAncestor].father);
        } else if (brute_tree_nodes[id][sAncestor].height < brute_tree_nodes[id][tAncestor].height) {
            if (static_cast<int>(sAncestor) == brute_tree_nodes[id][tAncestor].father)
                return make_pair(sAncestor, sAncestor);
            tAncestor = static_cast<unsigned int>(brute_tree_nodes[id][tAncestor].father);
        } else {
            if (brute_tree_nodes[id][sAncestor].father == -1 || sAncestor == tAncestor)
                return make_pair(sAncestor, tAncestor);
            sAncestor = static_cast<unsigned int>(brute_tree_nodes[id][sAncestor].father);
            tAncestor = static_cast<unsigned int>(brute_tree_nodes[id][tAncestor].father);
        }
    }
}

vector<pair<unsigned, double> > SegmentGraph::brute_down_top_search(unsigned down, unsigned up, unsigned id) {
    vector<pair<unsigned, double>> result;
    unordered_map<unsigned, double> distmap;
    distmap[down] = 0;
    auto anc = down;
    while (anc != up) {
        unordered_map<unsigned int, double> tmp;
        SegAttr &edges = brute_tree_nodes[id][anc].edges;
        double dist = FLOAT_IVF;
        if (distmap.find(anc) != distmap.end())
            dist = distmap[anc];
        for (auto &u: edges) {
            for (auto &medge: u.second.list) {
                if (medge.check_in_segment(1, 1)) {
                    tmp.insert(make_pair(u.first, medge.weight + dist));
                    break;
                }
            }
            if (distmap.find(u.first) != distmap.end()) {
                if (tmp.find(u.first) != distmap.end())
                    tmp[u.first] = min(tmp[u.first], distmap[u.first]);
                else tmp[u.first] = distmap[u.first];
            }
        }
        distmap = std::move(tmp);
        anc = static_cast<unsigned>(brute_tree_nodes[id][anc].father);
    }
    for (auto &p: distmap)
        result.emplace_back(p);
    return result;
}


double SegmentGraph::brute_range_decompose_query(unsigned u, unsigned v, unsigned left_range, unsigned right_range) {
    unsigned id = get_range_id(left_range, right_range);
    pair<unsigned, unsigned> lca = brute_qLCA(u, v, id);
    vector<pair<unsigned int, double>> s_result, t_result;
    s_result = brute_down_top_search(u, lca.first, id);
    t_result = brute_down_top_search(v, lca.second, id);
    double distance = FLOAT_IVF;
    for (auto &a: s_result) {
        if (a.second == FLOAT_IVF) continue;
        for (auto &b: t_result) {
            if (b.second == FLOAT_IVF) continue;
            if (a.first == b.first) distance = min(distance, a.second + b.second);
            else {
                distance = min(distance, a.second + b.second + brute_cp_distv1(a.first, b.first, id));
            }
        }
    }
    return distance;
}
