//
// Created by Mingyu on 23-10-9.
//
#include "Graph.h"
#include "utils.h"
#include "Segment.h"

#ifndef SEGMENTSTP_SEGMENTGRAPH_H

union PLLIndex {
    PLLIndex() {}

    ~PLLIndex() {
        std::vector<std::unordered_map<unsigned, double> >().swap(Label);
        std::vector<int>().swap(node_order_);
        std::vector<int>().swap(invert_node_order_);
    }

    std::vector<int> node_order_;
    std::vector<int> invert_node_order_;
    std::vector<std::unordered_map<unsigned, double> > Label;
};

union H2HIndex {
    H2HIndex() {}

    ~H2HIndex() {
        std::vector<unsigned>().swap(rank_);
        std::vector<TreeNode>().swap(H2H_Tree_);
        std::vector<unsigned>().swap(toRMQ);
        std::vector<std::vector<unsigned>>().swap(RMQIndex);
        std::vector<int>().swap(node_order_);
        std::vector<int>().swap(invert_node_order_);
    }

    std::vector<int> node_order_;
    std::vector<int> invert_node_order_;
    std::vector<unsigned> rank_;
    std::vector<TreeNode> H2H_Tree_;
    std::vector<unsigned> toRMQ;
    std::vector<std::vector<unsigned>> RMQIndex;
};

#define SEGMENTSTP_SEGMENTGRAPH_H


class SegmentGraph : public GraphIndex {
public:

    /*
     * file operation
     */

    void save_segraph_index(char *filename);

    void load_segraph_index(char *filename);

    void save_segraph_index_add(char *filename);

    void load_segraph_index_add(char *filename);

    void save_seg_pll_index(char *filename);

    void load_seg_pll_index(char *filename);

    /*
     * tree analysis
     */

    void calculate_height();

    /*
     * segment operations
     */

    unsigned segment_mapping(unsigned left_range, unsigned right_range);

    unsigned get_range_id(unsigned left_range, unsigned right_range);

    void range_filter_search(unsigned left_range, unsigned right_range);

    void random_walk_query(unsigned &u, unsigned &v, unsigned &left, unsigned &right, unsigned length) {
        unsigned L = 0, R = 0;
        u = rand() % node_num_;
        v = u;
        unsigned count = 1000*length;
        while (R - L + 1 <= length && count) {
            unsigned rd_ins = rand() % weight_graph_[v].size();
            if (L == 0) {
                L = attr_graph_[v][rd_ins];
                R = L;
            } else {
                L = std::min(L, attr_graph_[v][rd_ins]);
                R = std::max(R, attr_graph_[v][rd_ins]);
                left = L;
                right = R;
            }
            if (R - L + 1 <= length) v = weight_graph_[v][rd_ins].first;
            count--;
        }

    }


    /*
     * data modify
     */

    void add_random_range_attr(unsigned left_range, unsigned right_range);

    /*
     * range h2h index
     */

    void build_range_tree_decompose(unsigned type = 1);

    void build_range_tree_parallel();

    void build_range_tree_index_parallel();

    void build_range_tree_parallel0();

    void build_range_tree_index_parallel0();

    void BU_tree_index_shortest_path_task_para(unsigned id);

    void TD_tree_index_shortest_path_task_para(unsigned id);

    void range_tree_edge_decomposition();

    void range_tree_node_decomposition();

    void build_range_tree();

    void build_range_tree_index();

    inline const SegList &cp_distv2(unsigned s, unsigned t) {
        if (node_order_[s] < node_order_[t]) {
            swap(s, t);
        }
        return tree_nodes_[t].edges[s];
    }

    inline double cp_distv1(unsigned s, unsigned t, unsigned l, unsigned r) {
        assert(s != t);
        if (node_order_[s] < node_order_[t])
            swap(s, t);
        return tree_nodes_[t].edges[s].dist(l, r);
    }

    inline double brute_cp_distv1(unsigned s, unsigned t, unsigned id) {
        assert(s != t);
        if (preserve_order[id][s] < preserve_order[id][t])
            swap(s, t);
        return brute_tree_nodes[id][t].edges[s].dist(1, 1);
    }


    pair<unsigned int, unsigned int> qLCA(unsigned s, unsigned t);

    vector<pair<unsigned, double> >
    down_top_search(unsigned down, unsigned up, unsigned left_range, unsigned right_range);

    double range_tree_decompose_query(unsigned u, unsigned v, unsigned left_range, unsigned right_range);


    pair<unsigned int, unsigned int> brute_qLCA(unsigned s, unsigned t, unsigned id);

    vector<pair<unsigned, double> > brute_down_top_search(unsigned down, unsigned up, unsigned id);

    double brute_range_decompose_query(unsigned u, unsigned v, unsigned left_range, unsigned right_range);

    /*
     * order by tree width; tree height; seg list size; ave seg list size;
     */
    std::vector<double> tree_statistics_analysis();

    /*
     * range PLL index
     */

    void build_range_PLL_index(unsigned left_range, unsigned right_range);


    void build_brute_range_tree();

    void build_brute_range_tree_and_save(char *filename);


    /*
     * shortest path with dijkstra
     */

    double dijkstra_range_path(unsigned u, unsigned v, WeightGraph &cur_graph_, Segment seg);

    /*
     * shortest path with bi-dijkstra
     */

    double bi_dijkstra_range_path(unsigned u, unsigned v, WeightGraph &cur_graph_, Segment seg);

    /*
     * shortest path with PLL range index
     */

    double range_PLL_query(unsigned u, unsigned v, unsigned left_range, unsigned right_range);


    std::vector<std::vector<unsigned>> seg_id_map_;
    unsigned Left_Bound_, Right_Bound_;
    WeightGraph preserve_graph;
    AttrGraph attr_graph_, preserve_attr_;
    SegHop seg_hop_;
    std::vector<SegTreeNode> tree_nodes_, tmp_tree_nodes_;
    std::vector<unsigned> roots;
    std::vector<std::vector<int> > preserve_order;

    unordered_map<unsigned, vector<pair<unsigned, SegList>>> hopIndex;

    //for tmp results
    unordered_map<unsigned, SegList> tmp_dist;
    vector<unordered_map<unsigned int, vector<unsigned int>>> short_cut_dependency;
    vector<vector<vector<unsigned >>> t2h2v;
    vector<vector<pair<unsigned, unsigned>>> h2e;
    unsigned parallel_count = 0, height_threshold;
    mutex *lock;


    //for brute n2 range
    std::vector<std::vector<SegTreeNode>> brute_tree_nodes;

    /*
     * test data generate
     */

    void generate_small_test_graph(unsigned num, unsigned range) {
        node_num_ = num;
        weight_graph_.resize(num);
        attr_graph_.resize(num);
        srand(123);
        Left_Bound_ = 1;
        Right_Bound_ = range;
        std::vector<std::vector<bool>> dist_mp;
        dist_mp.resize(num);
        for (int i = 0; i < num; i++) dist_mp[i].resize(num);
        for (int i = 0; i < num; i++) {
            for (int j = 0; j < num; j++) {
                dist_mp[i][j] = false;
            }
        }

        for (int i = 1; i < num; i++) {
            unsigned next_node = rand() % (i + 1);
            double w = rand() % 10 + 1;
            if (dist_mp[i][next_node]) continue;
            dist_mp[i][next_node] = true;
            dist_mp[next_node][i] = true;
            weight_graph_[i].emplace_back(next_node, w);
            weight_graph_[next_node].emplace_back(i, w);
            attr_graph_[i].push_back(i % range + 1);
            attr_graph_[next_node].push_back(i % range + 1);
        }
        for (int i = 1; i < num; i++) {
            unsigned next_node = rand() % (i + 1);
            double w = rand() % 10 + 1;
            if (dist_mp[i][next_node]) continue;
            dist_mp[i][next_node] = true;
            dist_mp[next_node][i] = true;
            weight_graph_[i].emplace_back(next_node, w);
            weight_graph_[next_node].emplace_back(i, w);
            attr_graph_[i].push_back(i % range + 1);
            attr_graph_[next_node].push_back(i % range + 1);
        }


        for (int i = 0; i < num; i++) {
            for (int j = 0; j < weight_graph_[i].size(); j++) {
                std::cerr << i << " " << weight_graph_[i][j].first << " " << weight_graph_[i][j].second << endl;
            }
        }
    }

    void print_tree_node() {
        queue<unsigned> Q;
        for (auto root: roots) {
            std::cerr << "tree root:: " << root << endl;
            Q.push(root);
        }
        while (!Q.empty()) {
            unsigned u = Q.front();
            Q.pop();
            for (auto v: tree_nodes_[u].ch) {
                std::cerr << u << " " << v << endl;
                Q.push(v);
            }
        }
    }
};

#endif //SEGMENTSTP_SEGMENTGRAPH_H
