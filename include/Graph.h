//
// Created by Mingyu on 23-9-24.
//
#include "utils.h"

#ifndef SEGMENTSTP_GRAPH_H
typedef std::vector<std::vector<unsigned>> NormalGraph;
typedef std::vector<std::vector<std::pair<unsigned, double>>> WeightGraph;
typedef std::priority_queue<std::pair<double, unsigned>, std::vector<std::pair<double, unsigned>>, std::greater<> > NodeQueue;
typedef std::priority_queue<std::pair<unsigned, unsigned>, std::vector<std::pair<unsigned, unsigned>>, std::greater<> > PairQueue;
#define SEGMENTSTP_GRAPH_H

class GraphIndex {
public:
    /*
     * data load
     */

    void load_normal_graph(char *filename);

    void load_weight_graph(char *filename);

    void load_PLL_index(char *filename, const std::string &file_type);

    void save_PLL_index(char *filename, const std::string &file_type);

    void load_CH_index(char *filename, const std::string &file_type);

    void save_CH_index(char *filename, const std::string &file_type);

    void load_H2H_index(char *filename, const std::string &file_type);

    void save_H2H_index(char *filename, const std::string &file_type);

    void load_order_file(char *filename, const std::string &file_type);

    void save_order_file(char *filename, const std::string &file_type);

    /*
     *  contraction Hierarchies decomposition
     */

    void eliminate_E_order(unsigned u, unsigned v);

    void insert_E_order(unsigned u, unsigned v, double w);

    void CH_decomposition();

    /*
     * H2H construct
     */

    int match(std::vector<std::pair<unsigned, double>> &nodes);

    void build_tree();

    void build_tree_index();

    void build_H2H_index();

    void build_tree_index_dfs(unsigned u,std::vector<unsigned>& list);

    /*
     * PLL construct
     */

    void sort_node_by_degree();

    void dijkstra_prune(unsigned u, std::vector<std::pair<unsigned, double>> &vp);

    void PLL_construct();

    double PLL_SupNode_query(unsigned u, unsigned v, std::vector<unsigned> &SupNode, double &dist);

    /*
     * RMQ LCA
     */

    void build_RMQ();

    void build_RMQ_dfs_order(unsigned u);

    unsigned LCA_query(unsigned u,unsigned v);



    /*
     * shortest path without index
     */

    static int bfs_shortest_path(unsigned u, unsigned v, WeightGraph &cur_graph_);

    static double dijkstra_shortest_path(unsigned u, unsigned v, WeightGraph &cur_graph_);

    static double bi_dijkstra_shortest_path(unsigned u, unsigned v, WeightGraph &cur_graph_);


    /*
     * shortest path with CH
     */

    double CH_short_dist_query(unsigned u, unsigned v);

    /*
     * shortest path with H2H
     */

    double H2H_short_dist_query(unsigned u, unsigned v);

    /*
     * shortest path with PLL
     */

    double PLL_short_dist_query(unsigned u, unsigned v);


    /*
     * graph data
     */

    unsigned node_num_, edge_num_;
    WeightGraph weight_graph_;

    /*
     * node order
     */

    std::vector<int> node_order_;
    std::vector<int> invert_node_order_;

    /*
     * PLL data
     */
    std::vector<std::unordered_map<unsigned, double> > Label;

    /*
     * CH data
     */

    std::vector<std::unordered_map<unsigned, double> > E;
    std::vector<std::vector<std::pair<unsigned, double> > > CH_graph_;
    std::vector<std::unordered_map<unsigned, std::vector<unsigned> > > CH_short_cut_;

    /*
     * H2H data
     */

    unsigned height_max_;
    std::vector<std::vector<unsigned >> node_to_tree_map_; //one vertex exist in those tree nodes (nodeID--->tree node rank)
    std::vector<unsigned> EulerSeq; //prepare for the LCA calculation

    std::vector<unsigned> rank_;
    std::vector<TreeNode> H2H_Tree_;
    std::vector<unsigned> toRMQ;
    std::vector<std::vector<unsigned>> RMQIndex;

    /*
     * index status
     */
    bool is_load_pll = false, is_load_ch = false, is_load_h2h = false;

    unsigned num_threads = 80;

};


#endif //SEGMENTSTP_GRAPH_H
