//
// Created by mingyu on 23-10-11.
//

#include "SegmentGraph.h"
#include "Graph.h"
#include "utils.h"
using namespace std;

unsigned SegmentGraph::segment_mapping(unsigned left_range, unsigned right_range) {
    seg_id_map_.resize(right_range+1);
    unsigned cnt = 0;
    for (int i = 0; i <= right_range; i++) seg_id_map_[i].resize(right_range+1);
    for (unsigned i = left_range; i <=right_range ; i++) {
        for (unsigned j = i; j <= right_range; j++) {
            seg_id_map_[i][j] = cnt;
            cnt++;
        }
    }
    return cnt;
}


void SegmentGraph::calculate_height() {
    queue<unsigned> Q;
    for (auto v: roots) {
        tree_nodes_[v].height = 0;
        for (auto v1: tree_nodes_[v].ch) {
            Q.push(v1);
            tree_nodes_[v1].height = 1;
        }
    }

    while (!Q.empty()) {
        auto v = Q.front();
        Q.pop();
        for (auto a: tree_nodes_[v].ch) {
            Q.push(a);
            tree_nodes_[a].height = tree_nodes_[v].height + 1;
        }
    }
}


unsigned SegmentGraph::get_range_id(unsigned int left_range, unsigned int right_range) {
    return seg_id_map_[left_range][right_range];
}


void SegmentGraph::range_filter_search(unsigned left_range, unsigned right_range) {
    std::cerr<<"filter:: "<<left_range<<" "<<right_range<<endl;
    weight_graph_.resize(node_num_);
    attr_graph_.resize(node_num_);
    for (int i = 0; i < preserve_graph.size(); i++) {
        weight_graph_[i].resize(0);
        attr_graph_[i].resize(0);
        for (int j = 0; j < preserve_graph[i].size(); j++) {
            if (!check_in_segment(preserve_attr_[i][j], std::make_pair(left_range, right_range))) continue;
            weight_graph_[i].push_back(preserve_graph[i][j]);
            attr_graph_[i].push_back(1);
        }
    }
}

/*
 * minimum degree decomposition
 */

void SegmentGraph::range_tree_edge_decomposition() {
    node_order_.resize(node_num_, -1);
    invert_node_order_.resize(node_num_);
    seg_hop_.resize(node_num_);
    tree_nodes_.resize(node_num_);
    std::vector<int> old_degree(node_num_), new_degree(node_num_);

    auto cmp = [&old_degree](const unsigned &a, const unsigned &b) -> bool {
        return old_degree[a] == old_degree[b] ? a < b : old_degree[a] < old_degree[b];
    };

    set<unsigned, decltype(cmp)> node_degree_set(cmp);

    for (int i = 0; i < node_num_; i++) {
        old_degree[i] = new_degree[i] = (int) weight_graph_[i].size();
        node_degree_set.insert(i);
        for (int j = 0; j < weight_graph_[i].size(); j++) {
            unsigned v = weight_graph_[i][j].first;
            double weight = weight_graph_[i][j].second;
            seg_hop_[i][v].emplace_back_seg(weight, attr_graph_[i][j], attr_graph_[i][j]);
        }
    }
    int node_order = 0;
    while (!node_degree_set.empty()) {
        unsigned cur_node = *node_degree_set.begin();
        node_degree_set.erase(node_degree_set.begin());

        while (old_degree[cur_node] != new_degree[cur_node]) {
            old_degree[cur_node] = new_degree[cur_node];
            node_degree_set.insert(cur_node);

            cur_node = *node_degree_set.begin();
            node_degree_set.erase(node_degree_set.begin());
        }
        invert_node_order_[node_order] = (int) cur_node;
        node_order_[cur_node] = node_order++;
        auto &v = seg_hop_[cur_node];
        std::vector<unsigned> cur_neighbor;
        for (const auto &u: v) {
            if (node_order_[u.first] == -1) {
                cur_neighbor.push_back(u.first);
            }
        }
        if (node_order % (node_num_ / 100 + 1) == 0)
            std::cerr << "current order:: " << node_order << " size:: " << cur_neighbor.size() << endl;
        vector<int> neighbor_degree_increase_cnt(cur_neighbor.size(), 0);

        for (int i = 0; i < cur_neighbor.size(); i++) {
            auto ipair = make_pair(cur_neighbor[i], v[cur_neighbor[i]]);
            neighbor_degree_increase_cnt[i] -= (int) ipair.second.size();
            for (int j = i + 1; j < cur_neighbor.size(); j++) {
                auto jpair = make_pair(cur_neighbor[j], v[cur_neighbor[j]]);
                unsigned inner_size = seg_hop_[jpair.first][ipair.first].size();
                seg_hop_[jpair.first][ipair.first].combine(ipair.second + jpair.second);
                seg_hop_[ipair.first][jpair.first] = seg_hop_[jpair.first][ipair.first];

                inner_size = seg_hop_[jpair.first][ipair.first].size() - inner_size;
                neighbor_degree_increase_cnt[i] += (int) inner_size;
                neighbor_degree_increase_cnt[j] += (int) inner_size;
            }
            new_degree[ipair.first] += neighbor_degree_increase_cnt[i];
            if (neighbor_degree_increase_cnt[i] < 0) {
                node_degree_set.erase(ipair.first);
                old_degree[ipair.first] = new_degree[ipair.first];
                node_degree_set.insert(ipair.first);
            }
        }

        for (auto i: cur_neighbor) {
            tree_nodes_[cur_node].edges.insert(std::move(make_pair(i, v[i])));
        }
        if (tree_nodes_[cur_node].edges.empty()) {
            roots.push_back(cur_node);
        }
    }
}


void SegmentGraph::range_tree_node_decomposition() {
    vector<vector<unsigned>> degree2nodeQ;//degree->a list of vertex of this degree
    vector<pair<unsigned, unsigned >> vPosition(node_num_);//(degree,idx)

    for (unsigned v = 0; v < weight_graph_.size(); ++v) {
        unsigned degree = weight_graph_[v].size();
        if (degree >= degree2nodeQ.size()) {
            degree2nodeQ.resize(degree + 1);
        }
        vPosition[v] = make_pair(degree, degree2nodeQ[degree].size());
        degree2nodeQ[degree].push_back(v);
    }
    seg_hop_.resize(node_num_);
    for (int i = 0; i < node_num_; i++) {
        for (int j = 0; j < weight_graph_[i].size(); j++) {
            unsigned v = weight_graph_[i][j].first;
            double weight = weight_graph_[i][j].second;
            seg_hop_[i][v].emplace_back_seg(weight, attr_graph_[i][j], attr_graph_[i][j]);
        }
    }
    node_order_.resize(node_num_, -1);

    unsigned reduced_cnt = 0;
    unsigned int threshold_k = 0;
    vector<pair<unsigned int, unsigned int>> twidth_threshold2reduced_vcnt;
    tree_nodes_.resize(node_num_);
    unsigned mindegree = 0;
    int vorder = 0;

    while (true) {
        while (mindegree < degree2nodeQ.size() && degree2nodeQ[mindegree].empty()) {
            if (mindegree == threshold_k) {
                if (reduced_cnt != 0) {
                    twidth_threshold2reduced_vcnt.emplace_back(mindegree, reduced_cnt);
                    reduced_cnt = 0;
                }
                threshold_k++;
            }
            mindegree++;
        }

        if (degree2nodeQ.size() == mindegree)
            break;
        unsigned vid = degree2nodeQ[mindegree].back();
        degree2nodeQ[mindegree].pop_back();
        node_order_[vid] = vorder++;

        //adding short cuts
        auto &v = seg_hop_[vid];
        std::vector<unsigned> cur_neighbor;
        for (const auto &u: v) {
            if (node_order_[u.first] == -1) {
                cur_neighbor.push_back(u.first);
            }
        }
        if (vorder % (node_num_ / 100 + 1) == 0)
            std::cerr << "current order:: " << vorder << " size:: " << cur_neighbor.size() << endl;
        vector<int> neighbor_degree_increase_cnt(cur_neighbor.size(), -1);
        //add short cuts
        for (unsigned int i = 0; i < cur_neighbor.size(); ++i) {
            for (unsigned int j = i + 1; j < cur_neighbor.size(); ++j) {
                auto ipair = make_pair(cur_neighbor[i], v[cur_neighbor[i]]);
                auto jpair = make_pair(cur_neighbor[j], v[cur_neighbor[j]]);
                if (seg_hop_[ipair.first].find(jpair.first) == seg_hop_[ipair.first].end()) {
                    neighbor_degree_increase_cnt[i]++;
                    neighbor_degree_increase_cnt[j]++;
                }
                seg_hop_[jpair.first][ipair.first].combine(ipair.second + jpair.second);
                seg_hop_[ipair.first][jpair.first] = seg_hop_[jpair.first][ipair.first];
            }
        }
        //update vPosition
        for (unsigned int i = 0; i < cur_neighbor.size(); ++i) {

            if (neighbor_degree_increase_cnt[i] != 0) {
                unsigned int &x = cur_neighbor[i];
                pair<unsigned, unsigned> &p = vPosition[x];
                //swap and delete
                degree2nodeQ[p.first][p.second] = degree2nodeQ[p.first].back();
                vPosition[degree2nodeQ[p.first].back()].second = p.second;
                degree2nodeQ[p.first].pop_back();
                //place in a new position
                p.first += static_cast<unsigned int>(neighbor_degree_increase_cnt[i]);
                if (p.first >= degree2nodeQ.size()) {
                    degree2nodeQ.resize(p.first + 1);
                }
                mindegree = min(mindegree, p.first);
                p.second = degree2nodeQ[p.first].size();
                degree2nodeQ[p.first].push_back(x);
            }
        }
        //building tnodes
        for (auto i: cur_neighbor) {
            tree_nodes_[vid].edges.insert(std::move(make_pair(i, v[i])));
        }
        if (tree_nodes_[vid].edges.empty()) {
            roots.push_back(vid);
        }
    }
}


void SegmentGraph::build_range_tree() {
    for (unsigned i = 0; i < tree_nodes_.size(); i++) {
        if (i % (node_num_ / 100 + 1) == 0) std::cerr << "current range tree nodes:: " << i << endl;
        auto &node = tree_nodes_[i];
        auto &edges = node.edges;
        if (!edges.empty()) {
            unsigned next_node = edges.begin()->first;
            unsigned idx = next_node;
            int order = node_order_[next_node];
            for (auto &u: edges) {
                next_node = u.first;
                if (node_order_[next_node] != -1) {
                    if (order == -1 || node_order_[next_node] < order) {
                        idx = next_node;
                        order = node_order_[next_node];
                    }
                }
            }
            node.father = (int) idx;
            tree_nodes_[idx].ch.push_back(i);
        } else {
            node.father = -1;
        }
    }
}

/*
 * bfs build tree
 */
void SegmentGraph::build_range_tree_index() {
    unsigned cnt = 0;
    std::queue<unsigned> Q;
    for (auto u: roots) {
        tree_nodes_[u].height = 0;
        for (auto v: tree_nodes_[u].ch) {
            Q.push(v);
            tree_nodes_[v].height = 1;
        }
    }
    while (!Q.empty()) {
        unsigned v = Q.front();
        Q.pop();
        cnt++;
        if (cnt % (node_num_ / 100 + 1) == 0) std::cerr << "tree index node :: " << cnt << endl;
        for (auto a: tree_nodes_[v].ch) {
            Q.push(a);
            tree_nodes_[a].height = tree_nodes_[v].height + 1;
        }
        auto &edges = tree_nodes_[v].edges;

        for (auto neib1 = edges.begin(); neib1 != edges.end(); ++neib1) {
            for (auto neib2 = next(neib1); neib2 != edges.end(); ++neib2) {
                const SegList &d = cp_distv2(neib1->first, neib2->first);
                neib1->second.combine(neib2->second + d);
                neib2->second.combine(neib1->second + d);
            }
        }
    }
}

void SegmentGraph::build_range_tree_parallel() {

    vector<vector<unsigned>> degree2nodeQ;//degree->a list of vertex of this degree
    vector<pair<unsigned, unsigned >> vPosition(node_num_);//(degree,pos)

    for (unsigned v = 0; v < node_num_; ++v) {
        unsigned degree = weight_graph_[v].size();
        if (degree >= degree2nodeQ.size()) {
            degree2nodeQ.resize(degree + 1);
        }
        vPosition[v] = make_pair(degree, degree2nodeQ[degree].size());
        degree2nodeQ[degree].push_back(v);
    }

    short_cut_dependency.resize(node_num_);

    for (unsigned s = 0; s < weight_graph_.size(); ++s) {
        for (unsigned int i = 0; i < weight_graph_[s].size(); ++i) {
            auto &e = weight_graph_[s][i];
            short_cut_dependency[s][e.first].emplace_back(e.first);
        }
    }

    node_order_.resize(node_num_, -1);
    tree_nodes_.resize(node_num_);
    unsigned mindegree = 0;
    int vorder = 0;
    while (true) {
        while (mindegree < degree2nodeQ.size() && degree2nodeQ[mindegree].empty()) {
            mindegree++;
        }
        if (degree2nodeQ.size() == mindegree)
            break;
        unsigned vid = degree2nodeQ[mindegree].back();
        degree2nodeQ[mindegree].pop_back();
        node_order_[vid] = vorder++;
        //			mlog("reducting %u...", vid)
        //adding short cuts
        auto &v = short_cut_dependency[vid];
        std::vector<unsigned> cur_neighbor;
        for (const auto &u: v) {
            if (node_order_[u.first] == -1) {
                cur_neighbor.push_back(u.first);
            }
        }
        if (vorder % (node_num_ / 100 + 1) == 0)
            std::cerr << "current order:: " << vorder << " size:: " << cur_neighbor.size() << endl;
        vector<int> neighbor_degree_increase_cnt(cur_neighbor.size(), -1);
        //add short cuts
        for (unsigned int i = 0; i < cur_neighbor.size(); ++i) {
            for (unsigned int j = i + 1; j < cur_neighbor.size(); ++j) {
                auto ipair = make_pair(cur_neighbor[i], v[cur_neighbor[i]]);
                auto jpair = make_pair(cur_neighbor[j], v[cur_neighbor[j]]);
                if (short_cut_dependency[ipair.first].find(jpair.first) == short_cut_dependency[ipair.first].end()) {
                    neighbor_degree_increase_cnt[i]++;
                    neighbor_degree_increase_cnt[j]++;
                }
                short_cut_dependency[jpair.first][ipair.first].push_back(vid);
                short_cut_dependency[ipair.first][jpair.first].push_back(vid);
            }
        }

        //update vPosition
        for (unsigned int i = 0; i < cur_neighbor.size(); ++i) {

            if (neighbor_degree_increase_cnt[i] != 0) {
                unsigned int &x = cur_neighbor[i];
                pair<unsigned, unsigned> &p = vPosition[x];
                //swap and delete
                degree2nodeQ[p.first][p.second] = degree2nodeQ[p.first].back();
                vPosition[degree2nodeQ[p.first].back()].second = p.second;
                degree2nodeQ[p.first].pop_back();
                //place in a new position
                p.first += static_cast<unsigned int>(neighbor_degree_increase_cnt[i]);
                if (p.first >= degree2nodeQ.size()) {
                    degree2nodeQ.resize(p.first + 1);
                }
                mindegree = min(mindegree, p.first);
                p.second = degree2nodeQ[p.first].size();
                degree2nodeQ[p.first].push_back(x);
            }
        }

        //building tnodes
        for (auto i: cur_neighbor) {
            tree_nodes_[vid].edges[i];
        }

        if (tree_nodes_[vid].edges.empty()) {
            roots.push_back(vid);
        }
    }
    build_range_tree();
#pragma omp parallel num_threads(num_threads)
    {
#pragma omp single nowait
        {
            for (auto id: roots) {
#pragma omp task untied
                BU_tree_index_shortest_path_task_para(id);
            }
        }
    }
    short_cut_dependency.clear();
}

void SegmentGraph::BU_tree_index_shortest_path_task_para(unsigned id) {
    for (auto a: tree_nodes_[id].ch) {
#pragma omp task untied
        BU_tree_index_shortest_path_task_para(a);
    }
    /* first add its neighbors in original graph
     * then add the shortcuts
     */
    for (int j = 0; j < weight_graph_[id].size(); j++) {
        unsigned next = weight_graph_[id][j].first;
        unsigned attr = attr_graph_[id][j];
        double weight = weight_graph_[id][j].second;
        if (node_order_[next] > node_order_[id]) {
            tree_nodes_[id].edges[next].emplace_back_seg(weight, attr, attr);
        }
    }
#pragma omp taskwait
    int i = 0;
    for (auto &neibp: tree_nodes_[id].edges) {
        i = 0;
        if (short_cut_dependency[id][neibp.first].empty()) continue;
        if (short_cut_dependency[id][neibp.first][0] == neibp.first) i = 1;

        for (; i < short_cut_dependency[id][neibp.first].size(); ++i) {
            unsigned vid = short_cut_dependency[id][neibp.first][i];
            neibp.second.combine(tree_nodes_[vid].edges[id] + tree_nodes_[vid].edges[neibp.first]);
        }
    }
}

void SegmentGraph::build_range_tree_index_parallel() {
#pragma omp parallel num_threads(num_threads)
    {
#pragma omp for nowait
        for (unsigned int root: roots) {
            tree_nodes_[root].height = 0;
            for (auto v1: tree_nodes_[root].ch) {
                tree_nodes_[v1].height = 1;
#pragma omp task untied
                TD_tree_index_shortest_path_task_para(v1);
            }
        }
    }
}


void SegmentGraph::TD_tree_index_shortest_path_task_para(unsigned id) {
    auto &edges = tree_nodes_[id].edges;
    for (auto neib1 = edges.begin(); neib1 != edges.end(); ++neib1) {
        for (auto neib2 = next(neib1); neib2 != edges.end(); ++neib2) {
            unsigned s = neib1->first, t = neib2->first;
            if (node_order_[s] < node_order_[t]) {
                swap(s, t);
            }
            SegList &d = tmp_tree_nodes_[t].edges[s];
            neib1->second.combine(neib2->second + d);
            neib2->second.combine(neib1->second + d);
        }
    }
#pragma omp critical
    {
        parallel_count++;
        tmp_tree_nodes_[id] = tree_nodes_[id];
        if (parallel_count % (node_num_ / 100 + 1) == 0)
            std::cerr << "current parallel count:: " << parallel_count << endl;
    }
    for (auto a: tree_nodes_[id].ch) {
        tree_nodes_[a].height = tree_nodes_[id].height + 1;
#pragma omp task untied
        TD_tree_index_shortest_path_task_para(a);
    }
}


void SegmentGraph::build_range_tree_parallel0() {
    vector<vector<unsigned>> degree2nodeQ;//degree->a list of vertex of this degree
    vector<pair<unsigned, unsigned >> vPosition(node_num_);//(degree,pos)

    for (unsigned v = 0; v < node_num_; ++v) {
        unsigned degree = weight_graph_[v].size();
        if (degree >= degree2nodeQ.size()) {
            degree2nodeQ.resize(degree + 1);
        }
        vPosition[v] = make_pair(degree, degree2nodeQ[degree].size());
        degree2nodeQ[degree].push_back(v);
    }

    short_cut_dependency.resize(node_num_);

    for (unsigned s = 0; s < weight_graph_.size(); ++s) {
        for (unsigned int i = 0; i < weight_graph_[s].size(); ++i) {
            auto &e = weight_graph_[s][i];
            short_cut_dependency[s][e.first].emplace_back(e.first);
        }
    }

    node_order_.resize(node_num_, -1);
    tree_nodes_.resize(node_num_);
    unsigned mindegree = 0;
    int vorder = 0;
    while (true) {
        while (mindegree < degree2nodeQ.size() && degree2nodeQ[mindegree].empty()) {
            mindegree++;
        }
        if (degree2nodeQ.size() == mindegree)
            break;
        unsigned vid = degree2nodeQ[mindegree].back();
        degree2nodeQ[mindegree].pop_back();
        node_order_[vid] = vorder++;
        //			mlog("reducting %u...", vid)
        //adding short cuts
        auto &v = short_cut_dependency[vid];
        std::vector<unsigned> cur_neighbor;

        for (const auto &u: v) {
            if (node_order_[u.first] == -1) {
                cur_neighbor.push_back(u.first);
            }
        }
        if (vorder % (node_num_ / 100 + 1) == 0)
            std::cerr << "current order:: " << vorder << " size:: " << cur_neighbor.size() << endl;
        vector<int> neighbor_degree_increase_cnt(cur_neighbor.size(), -1);
        //add short cuts
        for (unsigned int i = 0; i < cur_neighbor.size(); ++i) {
            for (unsigned int j = i + 1; j < cur_neighbor.size(); ++j) {
                auto ipair = make_pair(cur_neighbor[i], v[cur_neighbor[i]]);
                auto jpair = make_pair(cur_neighbor[j], v[cur_neighbor[j]]);
                if (short_cut_dependency[ipair.first].find(jpair.first) == short_cut_dependency[ipair.first].end()) {
                    neighbor_degree_increase_cnt[i]++;
                    neighbor_degree_increase_cnt[j]++;
                }
                short_cut_dependency[jpair.first][ipair.first].push_back(vid);
                short_cut_dependency[ipair.first][jpair.first].push_back(vid);
            }
        }

        //update vPosition
        for (unsigned int i = 0; i < cur_neighbor.size(); ++i) {
            if (neighbor_degree_increase_cnt[i] != 0) {
                unsigned int &x = cur_neighbor[i];
                pair<unsigned, unsigned> &p = vPosition[x];
                //swap and delete
                degree2nodeQ[p.first][p.second] = degree2nodeQ[p.first].back();
                vPosition[degree2nodeQ[p.first].back()].second = p.second;
                degree2nodeQ[p.first].pop_back();
                //place in a new position
                p.first += static_cast<unsigned int>(neighbor_degree_increase_cnt[i]);
                if (p.first >= degree2nodeQ.size()) {
                    degree2nodeQ.resize(p.first + 1);
                }
                mindegree = min(mindegree, p.first);
                p.second = degree2nodeQ[p.first].size();
                degree2nodeQ[p.first].push_back(x);
            }
        }

        //building tnodes
        for (auto i: cur_neighbor) {
            tree_nodes_[vid].edges[i];
        }

        if (tree_nodes_[vid].edges.empty()) {
            roots.push_back(vid);
        }
    }
    build_range_tree();
    calculate_height();
    std::cerr << "parallel begin" << endl;
    vector<unsigned> vertex_height(weight_graph_.size(), 0);
    std::function<unsigned(unsigned)> cal_height = [&vertex_height, &cal_height, this](unsigned i) -> unsigned {
        for (auto c: tree_nodes_[i].ch) {
            vertex_height[i] = max(vertex_height[i], cal_height(c));
        }
        return vertex_height[i] + 1;
    };
    unsigned max_height = 0;
    for (auto cvid: roots) max_height = max(cal_height(cvid), max_height);

    h2e.resize(max_height + 1);
#pragma omp parallel num_threads(num_threads)
    {
        vector<vector<unsigned >> h2v(max_height + 1);
        //create h2v mapping
#pragma omp for nowait
        for (unsigned vid = 0; vid < vertex_height.size(); vid++) {
            h2v[vertex_height[vid]].push_back(vid);
        }

        //add shortcuts from original graph
#pragma omp for
        for (unsigned vid = 0; vid < weight_graph_.size(); vid++) {
            for (unsigned j = 0; j < weight_graph_[vid].size(); j++) {
                unsigned attr = attr_graph_[vid][j];
                auto ve = weight_graph_[vid][j];
                if (node_order_[ve.first] > node_order_[vid]) {//only add shortcuts from higher order
                    tree_nodes_[vid].edges[ve.first].emplace_back_seg(ve.second, attr, attr);
                }
            }
        }
        //first phase. process large amount of vertexes with small degrees and low height
        for (unsigned int height = 1; height <= min(height_threshold, max_height); ++height) {
            for (unsigned vid: h2v[height]) {
                //add from contributor
                for (auto &neibp: tree_nodes_[vid].edges) {
                    for (int i = 0; i < short_cut_dependency[vid][neibp.first].size(); ++i) {
                        unsigned id = short_cut_dependency[vid][neibp.first][i];
                        LockGuard guard(lock[id]);
                        neibp.second.combine(tree_nodes_[id].edges[vid] + tree_nodes_[id].edges[neibp.first]);
                    }
                }
            }
#pragma omp barrier
        }

        //second phase for small number of vertexes with high degree and large height
        //2.1 move local data out
#pragma omp critical(c1)
        t2h2v.emplace_back(std::move(h2v));
#pragma omp barrier
        //2.2 unfold edges to be calculated
#pragma omp for schedule(static, 4)//time this code and try other schedule type
        for (unsigned height = height_threshold + 1; height <= max_height; height++) {
            vector<pair<unsigned, unsigned >> es;
            for (auto &_h2v: t2h2v) {
                for (auto &_vid: _h2v[height]) {
                    for (auto &neibp: tree_nodes_[_vid].edges) es.emplace_back(_vid, neibp.first);
                }
            }
            h2e[height] = std::move(es);
        }
        //2.3 from params.height_threashold to max_height, parallel edges in each height
#pragma omp for
        for (unsigned int height = height_threshold + 1; height <= max_height; ++height) {
            auto &es = h2e[height];
            for (int i = 0; i < es.size(); i++) {
                auto &vid = es[i].first;
                auto &vidt = es[i].second;
                SegList &attr = tree_nodes_[es[i].first].edges[h2e[height][i].second];
                for (int j = 0; j < short_cut_dependency[vid][vidt].size(); ++j) {
                    unsigned id = short_cut_dependency[vid][vidt][j];
                    LockGuard guard(lock[id]);
                    attr.combine(tree_nodes_[id].edges[vid] + tree_nodes_[id].edges[vidt]);
                }
            }
        }
        //2.4 done
    }
    short_cut_dependency.clear();
}


void SegmentGraph::build_range_tree_index_parallel0() {
    unsigned cnt = 0;
    std::queue<unsigned> Q;
    for (auto u: roots) {
        tree_nodes_[u].height = 0;
        for (auto v: tree_nodes_[u].ch) {
            Q.push(v);
            tree_nodes_[v].height = 1;
        }
    }
    while (!Q.empty()) {
        unsigned v = Q.front();
        Q.pop();
        cnt++;
        if (cnt % (node_num_ / 100 + 1) == 0) std::cerr << "tree index node :: " << cnt << endl;
        for (auto a: tree_nodes_[v].ch) {
            Q.push(a);
            tree_nodes_[a].height = tree_nodes_[v].height + 1;
        }
        auto &edges = tree_nodes_[v].edges;

        if (edges.size() < 100) {
            for (auto neib1 = edges.begin(); neib1 != edges.end(); ++neib1) {
                for (auto neib2 = next(neib1); neib2 != edges.end(); ++neib2) {
                    const SegList &d = cp_distv2(neib1->first, neib2->first);
                    neib1->second.combine(neib2->second + d);
                    neib2->second.combine(neib1->second + d);
                }
            }
        } else {
            unsigned size = edges.size();
            SegAttr tmp_edges = tree_nodes_[v].edges;
            std::vector<SegAttr::iterator> tmp(size), task(size);
            unsigned iter_cnt = 0;
            for (auto neib1 = tmp_edges.begin(); neib1 != tmp_edges.end(); ++neib1) {
                tmp[iter_cnt] = neib1;
                iter_cnt++;
            }
            iter_cnt = 0;
            for (auto neib1 = edges.begin(); neib1 != edges.end(); ++neib1) {
                task[iter_cnt] = neib1;
                iter_cnt++;
            }

#pragma omp parallel for num_threads(num_threads)
            for (int i = 0; i < size; i++) {
                for (int j = 0; j < size; j++) {
                    if (i == j) continue;
                    const SegList &d = cp_distv2(tmp[i]->first, tmp[j]->first);
                    task[i]->second.combine(tmp[j]->second + d);
                }
            }
        }
    }
}


void SegmentGraph::build_brute_range_tree() {
    unsigned length = segment_mapping(Left_Bound_, Right_Bound_);
    std::cerr<<"brute length:: "<<length<<std::endl;
    brute_tree_nodes.resize(length);
    preserve_order.resize(length);
    for (unsigned i = Left_Bound_; i <= Right_Bound_; i++) {
        for (unsigned j = i; j <= Right_Bound_; j++) {
            tree_nodes_.clear();
            node_order_.clear();
            seg_hop_.clear();
            roots.clear();
            unsigned range_id = get_range_id(i,j);
            range_filter_search(i,j);
            range_tree_node_decomposition();
            build_range_tree();
            build_range_tree_index();
            /************ DEBUGGING*/
            brute_tree_nodes[range_id] = tree_nodes_;
            preserve_order[range_id] = node_order_;
             /***********/
        }
    }
}



void SegmentGraph::build_brute_range_tree_and_save(char *filename) {
    std::ofstream fout(filename, std::ios::binary);
    unsigned length = segment_mapping(Left_Bound_, Right_Bound_);
    std::cerr<<"brute length:: "<<length<<std::endl;
    preserve_graph = weight_graph_;
    preserve_attr_ = attr_graph_;
    fout.write((char *) &length, sizeof(unsigned));
    for (unsigned i = Left_Bound_; i <= Right_Bound_; i++) {
        for (unsigned j = i; j <= Right_Bound_; j++) {
            tree_nodes_.clear();
            node_order_.clear();
            seg_hop_.clear();
            roots.clear();
            unsigned range_id = get_range_id(i,j);
            range_filter_search(i,j);
            range_tree_node_decomposition();
            build_range_tree();
            build_range_tree_index();
            for (auto &node: tree_nodes_) {
                fout.write((char *) &node.height, sizeof(int));
                fout.write((char *) &node.father, sizeof(int));
                unsigned cnt = node.ch.size();
                fout.write((char *) &cnt, sizeof(unsigned));
                fout.write((char *) node.ch.data(), sizeof(unsigned) * cnt);

                unsigned edge_cnt = node.edges.size();
                fout.write((char *) &edge_cnt, sizeof(unsigned));
                for (auto &u: node.edges) {
                    fout.write((char *) &u.first, sizeof(int));
                    unsigned list_cnt = u.second.list.size();
                    fout.write((char *) &list_cnt, sizeof(unsigned));
                    for (auto v: u.second.list) {
                        fout.write((char *) &v.weight, sizeof(double));
                        fout.write((char *) &v.L, sizeof(unsigned));
                        fout.write((char *) &v.R, sizeof(unsigned));
                    }
                }
            }
            fout.write((char *) node_order_.data(), sizeof(int) * node_num_);
        }
    }
    weight_graph_ = preserve_graph;
    attr_graph_ = preserve_attr_;
}


void SegmentGraph::build_range_tree_decompose(unsigned type) {
    if (type == 1) {
        range_tree_edge_decomposition();
        build_range_tree();
        build_range_tree_index();
    } else if (type == 2) {
        range_tree_node_decomposition();
        build_range_tree();
        build_range_tree_index();
    } else if (type == 3) {
        build_range_tree_parallel();
        std::cerr << "build range tree index::" << endl;
        tmp_tree_nodes_ = tree_nodes_;
        build_range_tree_index_parallel();
    } else if (type == 4) {
        build_range_tree_parallel0();
        std::cerr << "build range tree index::" << endl;
        build_range_tree_index_parallel0();
    } else if (type == 5) {
        preserve_graph = weight_graph_;
        preserve_attr_ = attr_graph_;
        build_brute_range_tree();
        weight_graph_ = preserve_graph;
        attr_graph_ = preserve_attr_;
    }
}




