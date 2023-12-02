//
// Created by Mingyu on 23-9-30.
//

#include "utils.h"
#include "Graph.h"
#include "SegmentGraph.h"

using namespace std;


void GraphIndex::load_normal_graph(char *filename) {
    std::ifstream fin(filename);
    unsigned nodes, edges, u, v;
    fin >> nodes >> edges;
    node_num_ = nodes + 1;
    edge_num_ = edges;
    weight_graph_.resize(nodes);
    for (int i = 0; i < edges; i++) {
        fin >> u >> v;
        weight_graph_[u].emplace_back(v, 1);
    }
    std::cerr << "node count:: " << nodes << " edge count:: " << edges << std::endl;
}


void GraphIndex::load_weight_graph(char *filename) {
    std::ifstream fin(filename);
    unsigned nodes, edges, u, v;
    double w;
    fin >> nodes >> edges;
    node_num_ = nodes + 1;
    edge_num_ = edges;
    weight_graph_.resize(node_num_);
    for (int i = 0; i < edges; i++) {
        fin >> u >> v >> w;
        weight_graph_[u].emplace_back(v, w);
    }
    std::cerr << "node count:: " << nodes << " edge count:: " << edges << std::endl;
}

void GraphIndex::save_PLL_index(char *filename, const std::string &file_type) {
    if (file_type == "bin") {
        std::ofstream out(filename, std::ios::binary);
        out.write((char *) &node_num_, sizeof(unsigned));
        for (int i = 0; i < node_num_; i++) {
            unsigned size = Label[i].size();
            out.write((char *) &i, sizeof(unsigned));
            out.write((char *) &node_order_[i], sizeof(int));
            out.write((char *) &size, sizeof(unsigned));
            for (auto u: Label[i]) {
                out.write((char *) &u.first, sizeof(unsigned));
                out.write((char *) &u.second, sizeof(double));
            }
        }
    } else if (file_type == "str") {
        std::ofstream out(filename);
        out.setf(ios::fixed, ios::floatfield);
        out.precision(4);
        out << node_num_ << endl;
        for (int i = 0; i < node_num_; i++) {
            out << i << " " << node_order_[i] << " " << Label[i].size() << endl;
            for (auto u: Label[i]) {
                out << " " << u.first << " " << u.second;
            }
            out << endl;
        }
    }
}

void GraphIndex::load_PLL_index(char *filename, const std::string &file_type) {
    if (file_type == "bin") {
        std::ifstream in(filename, std::ios::binary);
        in.read((char *) &node_num_, sizeof(unsigned));
        node_order_.resize(node_num_);
        Label.resize(node_num_);
        for (int i = 0; i < node_num_; i++) {
            unsigned id, size, hub;
            int order;
            double hub_dist;
            in.read((char *) &id, sizeof(unsigned));
            in.read((char *) &order, sizeof(int));
            in.read((char *) &size, sizeof(unsigned));
            node_order_[id] = order;
            for (int j = 0; j < size; j++) {
                in.read((char *) &hub, sizeof(unsigned));
                in.read((char *) &hub_dist, sizeof(double));
                Label[id].emplace(hub, hub_dist);
            }
        }
    } else if (file_type == "str") {
        std::ifstream in(filename);
        in >> node_num_;
        node_order_.resize(node_num_);
        Label.resize(node_num_);
        for (int i = 0; i < node_num_; i++) {
            unsigned id, size, hub;
            int order;
            double hub_dist;
            in >> id >> order >> size;
            node_order_[id] = order;
            for (int j = 0; j < size; j++) {
                in >> hub >> hub_dist;
                Label[id].emplace(hub, hub_dist);
            }
        }
    }
}


void GraphIndex::load_order_file(char *filename, const std::string &file_type) {
    if (file_type == "bin") {
        std::ifstream in(filename, ios::binary);
        in.read((char *) &node_num_, sizeof(unsigned));
        node_order_.resize(node_num_);
        invert_node_order_.resize(node_num_);
        int node_id, node_order;
        for (int i = 0; i < node_num_; i++) {
            in.read((char *) &node_id, sizeof(unsigned));
            in.read((char *) &node_order, sizeof(unsigned));
            node_order_[node_id] = node_order;
            if (node_order != -1) {
                invert_node_order_[node_order] = node_id;
            }
        }
    } else if (file_type == "str") {
        std::ifstream in(filename);
        in >> node_num_;
        node_order_.resize(node_num_);
        invert_node_order_.resize(node_num_);
        int node_id, node_order;
        for (int i = 0; i < node_num_; i++) {
            in >> node_id >> node_order;
            node_order_[node_id] = node_order;
            if (node_order != -1) {
                invert_node_order_[node_order] = node_id;
            }
        }
    }
}


void GraphIndex::save_order_file(char *filename, const std::string &file_type) {
    if (file_type == "bin") {
        std::ofstream out(filename, ios::binary);
        out.write((char *) &node_num_, sizeof(unsigned));
        for (int i = 0; i < node_num_; i++) {
            out.write((char *) &i, sizeof(unsigned));
            out.write((char *) &node_order_[i], sizeof(unsigned));
        }
    } else if (file_type == "str") {
        std::ofstream out(filename);
        out << node_num_ << endl;
        for (int i = 0; i < node_num_; i++) {
            out << i << " " << node_order_[i] << endl;
        }
    }
}


void GraphIndex::load_CH_index(char *filename, const std::string &file_type) {
    if (file_type == "bin") {
        std::ifstream in(filename, std::ios::binary);
        in.read((char *) &node_num_, sizeof(unsigned));
        node_order_.resize(node_num_);
        invert_node_order_.resize(node_num_);
        CH_graph_.resize(node_num_);
        for (int i = 0; i < node_num_; i++) {
            unsigned id, size, edge;
            int order, invert_order;
            double edge_weight;
            in.read((char *) &id, sizeof(unsigned));
            in.read((char *) &order, sizeof(int));
            in.read((char *) &invert_order, sizeof(int));
            in.read((char *) &size, sizeof(unsigned));
            node_order_[id] = order;
            invert_node_order_[id] = invert_order;
            for (int j = 0; j < size; j++) {
                in.read((char *) &edge, sizeof(unsigned));
                in.read((char *) &edge_weight, sizeof(double));
                CH_graph_[id].emplace_back(edge, edge_weight);
            }
        }
    } else if (file_type == "str") {
        std::ifstream in(filename);
        in >> node_num_;
        node_order_.resize(node_num_);
        invert_node_order_.resize(node_num_);
        CH_graph_.resize(node_num_);
        for (int i = 0; i < node_num_; i++) {
            unsigned id, size, edge;
            int order, invert_order;
            double edge_weight;
            in >> id >> order >> invert_order >> size;
            node_order_[id] = order;
            invert_node_order_[id] = invert_order;
            for (int j = 0; j < size; j++) {
                in >> edge >> edge_weight;
                CH_graph_[id].emplace_back(edge, edge_weight);
            }
        }
    }
}


void GraphIndex::save_CH_index(char *filename, const std::string &file_type) {
    if (file_type == "bin") {
        std::ofstream out(filename, std::ios::binary);
        out.write((char *) &node_num_, sizeof(unsigned));
        for (int i = 0; i < node_num_; i++) {
            unsigned size = CH_graph_[i].size();
            out.write((char *) &i, sizeof(unsigned));
            out.write((char *) &node_order_[i], sizeof(int));
            out.write((char *) &invert_node_order_[i], sizeof(int));
            out.write((char *) &size, sizeof(unsigned));
            for (auto u: CH_graph_[i]) {
                out.write((char *) &u.first, sizeof(unsigned));
                out.write((char *) &u.second, sizeof(double));
            }
        }
    } else if (file_type == "str") {
        std::ofstream out(filename);
        out.setf(ios::fixed, ios::floatfield);
        out.precision(4);
        out << node_num_ << endl;
        for (int i = 0; i < node_num_; i++) {
            out << i << " " << node_order_[i] << " " << invert_node_order_[i] << " " << CH_graph_[i].size() << endl;
            for (auto u: CH_graph_[i]) {
                out << " " << u.first << " " << u.second;
            }
            out << endl;
        }
    }
}


void GraphIndex::load_H2H_index(char *filename, const std::string &file_type) {
    if (file_type == "bin") {
        std::ifstream in(filename, std::ios::binary);
        in.read((char *) &node_num_, sizeof(unsigned));
        toRMQ.resize(node_num_);
        node_order_.resize(node_num_);
        invert_node_order_.resize(node_num_);
        rank_.resize(node_num_);
        for (int i = 0; i < node_num_; i++) {
            in.read((char *) &toRMQ[i], sizeof(unsigned));
            in.read((char *) &node_order_[i], sizeof(int));
            in.read((char *) &invert_node_order_[i], sizeof(int));
            in.read((char *) &rank_[i], sizeof(unsigned));
        }
        unsigned RMQ_size, index_size;
        in.read((char *) &RMQ_size, sizeof(unsigned));
        RMQIndex.resize(RMQ_size);
        for (int i = 0; i < RMQ_size; i++) {
            in.read((char *) &index_size, sizeof(unsigned));
            RMQIndex[i].resize(index_size);
            for (int j = 0; j < index_size; j++) {
                in.read((char *) &RMQIndex[i][j], sizeof(unsigned));
            }
        }
        H2H_Tree_.resize(node_num_);
        for (int i = 0; i < node_num_; i++) {
            unsigned tree_size;
            in.read((char *) &tree_size, sizeof(unsigned));
            H2H_Tree_[i].nodes_.resize(tree_size);
            for (int j = 0; j < tree_size; j++) {
                unsigned id;
                double dist;
                in.read((char *) &id, sizeof(unsigned));
                in.read((char *) &dist, sizeof(double));
                H2H_Tree_[i].nodes_[j] = std::make_pair(id, dist);
            }
            in.read((char *) &H2H_Tree_[i].unique_node, sizeof(unsigned));
            in.read((char *) &H2H_Tree_[i].height, sizeof(unsigned));
            in.read((char *) &H2H_Tree_[i].fa, sizeof(unsigned));
            in.read((char *) &H2H_Tree_[i].depth, sizeof(unsigned));
            unsigned pos_size;
            in.read((char *) &pos_size, sizeof(unsigned));
            H2H_Tree_[i].pos.resize(pos_size);
            for (int j = 0; j < pos_size; j++) {
                in.read((char *) &H2H_Tree_[i].pos[j], sizeof(unsigned));
            }
            unsigned dis_size;
            in.read((char *) &dis_size, sizeof(unsigned));
            H2H_Tree_[i].dis.resize(dis_size);
            for (int j = 0; j < dis_size; j++) {
                in.read((char *) &H2H_Tree_[i].dis[j], sizeof(double));
            }
            unsigned ch_size;
            in.read((char *) &ch_size, sizeof(unsigned));
            H2H_Tree_[i].ch.resize(ch_size);
            for (int j = 0; j < ch_size; j++) {
                in.read((char *) &H2H_Tree_[i].ch[j], sizeof(unsigned));
            }
        }
    } else if (file_type == "str") {
        std::ifstream in(filename);
        in >> node_num_;
        for (int i = 0; i < node_num_; i++) {
            in >> toRMQ[i] >> node_order_[i] >> invert_node_order_[i] >> rank_[i];
        }
        unsigned RMQ_size, index_size;
        in >> RMQ_size;
        RMQIndex.resize(RMQ_size);
        for (int i = 0; i < RMQ_size; i++) {
            in >> index_size;
            RMQIndex[i].resize(index_size);
            for (int j = 0; j < index_size; j++) {
                in >> RMQIndex[i][j];
            }
        }
        H2H_Tree_.resize(node_num_);
        for (int i = 0; i < node_num_; i++) {
            unsigned nodes_size, pos_size, dis_size, ch_size;
            in >> nodes_size;
            H2H_Tree_[i].nodes_.resize(nodes_size);
            for (int j = 0; j < nodes_size; j++) {
                in >> H2H_Tree_[i].nodes_[j].first >> H2H_Tree_[i].nodes_[j].second;
            }
            in >> H2H_Tree_[i].unique_node >> H2H_Tree_[i].height >> H2H_Tree_[i].fa >> H2H_Tree_[i].depth;
            in >> pos_size;
            H2H_Tree_[i].pos.resize(pos_size);
            for (int j = 0; j < pos_size; j++) in >> H2H_Tree_[i].pos[j];
            in >> dis_size;
            H2H_Tree_[i].pos.resize(dis_size);
            for (int j = 0; j < dis_size; j++) in >> H2H_Tree_[i].dis[j];
            in >> ch_size;
            H2H_Tree_[i].pos.resize(ch_size);
            for (int j = 0; j < ch_size; j++) in >> H2H_Tree_[i].ch[j];
        }
    }
}


void GraphIndex::save_H2H_index(char *filename, const std::string &file_type) {
    if (file_type == "bin") {
        std::ofstream out(filename, std::ios::binary);
        out.write((char *) &node_num_, sizeof(unsigned));
        for (int i = 0; i < node_num_; i++) {
            out.write((char *) &toRMQ[i], sizeof(unsigned));
            out.write((char *) &node_order_[i], sizeof(int));
            out.write((char *) &invert_node_order_[i], sizeof(int));
            out.write((char *) &rank_[i], sizeof(unsigned));
        }
        unsigned RMQ_size = RMQIndex.size(), index_size;
        out.write((char *) &RMQ_size, sizeof(unsigned));
        for (int i = 0; i < RMQ_size; i++) {
            index_size = RMQIndex[i].size();
            out.write((char *) &index_size, sizeof(unsigned));
            for (int j = 0; j < index_size; j++) {
                out.write((char *) &RMQIndex[i][j], sizeof(unsigned));
            }
        }
        for (int i = 0; i < node_num_; i++) {
            unsigned tree_size = H2H_Tree_[i].nodes_.size();
            out.write((char *) &tree_size, sizeof(unsigned));
            for (int j = 0; j < tree_size; j++) {
                unsigned id = H2H_Tree_[i].nodes_[j].first;
                double dist = H2H_Tree_[i].nodes_[j].second;
                out.write((char *) &id, sizeof(unsigned));
                out.write((char *) &dist, sizeof(double));
            }
            out.write((char *) &H2H_Tree_[i].unique_node, sizeof(unsigned));
            out.write((char *) &H2H_Tree_[i].height, sizeof(unsigned));
            out.write((char *) &H2H_Tree_[i].fa, sizeof(unsigned));
            out.write((char *) &H2H_Tree_[i].depth, sizeof(unsigned));
            unsigned pos_size = H2H_Tree_[i].pos.size();
            out.write((char *) &pos_size, sizeof(unsigned));
            for (int j = 0; j < pos_size; j++) {
                out.write((char *) &H2H_Tree_[i].pos[j], sizeof(unsigned));
            }
            unsigned dis_size = H2H_Tree_[i].dis.size();
            out.write((char *) &dis_size, sizeof(unsigned));
            for (int j = 0; j < dis_size; j++) {
                out.write((char *) &H2H_Tree_[i].dis[j], sizeof(double));
            }
            unsigned ch_size = H2H_Tree_[i].ch.size();
            out.write((char *) &ch_size, sizeof(unsigned));
            for (int j = 0; j < ch_size; j++) {
                out.write((char *) &H2H_Tree_[i].ch[j], sizeof(unsigned));
            }
        }
    } else if (file_type == "str") {
        std::ofstream out(filename);
        out.setf(ios::fixed, ios::floatfield);
        out.precision(4);
        out << node_num_ << endl;
        for (int i = 0; i < node_num_; i++) {
            out << toRMQ[i] << " " << node_order_[i] << " " << invert_node_order_[i] << " " << rank_[i] << endl;
        }
        out << RMQIndex.size() << endl;
        for (auto &i: RMQIndex) {
            out << i.size() << endl;
            for (unsigned int j: i) {
                out << j << " ";
            }
            out << endl;
        }
        for (int i = 0; i < node_num_; i++) {
            out << H2H_Tree_[i].nodes_.size() << endl;
            for (auto u: H2H_Tree_[i].nodes_) {
                out << u.first << " " << u.second << " ";
            }
            out << endl;
            out << H2H_Tree_[i].unique_node << endl;
            out << H2H_Tree_[i].height << endl;
            out << H2H_Tree_[i].fa << endl;
            out << H2H_Tree_[i].depth << endl;
            out << H2H_Tree_[i].pos.size() << endl;
            for (auto v: H2H_Tree_[i].pos) {
                out << v << " ";
            }
            out << endl;
            out << H2H_Tree_[i].dis.size() << endl;
            for (auto v: H2H_Tree_[i].dis) {
                out << v << " ";
            }
            out << endl;
            out << H2H_Tree_[i].ch.size() << endl;
            for (auto v: H2H_Tree_[i].ch) {
                out << v << " ";
            }
            out << endl;
        }
    }
}


void SegmentGraph::save_segraph_index(char *filename) {
    std::ofstream fout(filename, std::ios::binary);
    unsigned check = 0;
    for (auto &node: tree_nodes_) {
        check++;
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


void SegmentGraph::load_segraph_index(char *filename) {
    std::ifstream fin(filename, std::ios::binary);
    tree_nodes_.resize(node_num_);
    node_order_.resize(node_num_);
    unsigned cnt;
    unsigned check = 0;
    for (auto &node: tree_nodes_) {
        check++;
        fin.read((char *) &node.height, sizeof(int));
        fin.read((char *) &node.father, sizeof(int));
        fin.read((char *) &cnt, sizeof(unsigned));
        node.ch.resize(cnt);
        fin.read((char *) node.ch.data(), sizeof(unsigned) * cnt);
        unsigned edge_cnt;
        fin.read((char *) &edge_cnt, sizeof(unsigned));
        for (int i = 0; i < edge_cnt; i++) {
            int id;
            fin.read((char *) &id, sizeof(int));
            unsigned list_cnt;
            fin.read((char *) &list_cnt, sizeof(unsigned));
            for (int j = 0; j < list_cnt; j++) {
                double weight;
                unsigned L, R;
                fin.read((char *) &weight, sizeof(double));
                fin.read((char *) &L, sizeof(unsigned));
                fin.read((char *) &R, sizeof(unsigned));
                node.edges[id].emplace_back_seg(weight, L, R);
            }
        }
    }
    fin.read((char *) node_order_.data(), sizeof(int) * node_num_);
}

void SegmentGraph::load_seg_pll_index(char *filename) {
    std::ifstream fin(filename, std::ios::binary);
    for (int i = 0; i < node_num_; i++) {
        unsigned size;
        fin.read((char *) &size, sizeof(unsigned));
        for (int j = 0; j < size; j++) {
            unsigned id, list_size;
            fin.read((char *) &id, sizeof(unsigned));
            fin.read((char *) &list_size, sizeof(unsigned));
            hopIndex[i].emplace_back(id, SegList());
            for (int k = 0; k < list_size; k++) {
                double weight;
                unsigned L, R;
                fin.read((char *) &weight, sizeof(double));
                fin.read((char *) &L, sizeof(unsigned));
                fin.read((char *) &R, sizeof(unsigned));
                hopIndex[i][j].second.emplace_back_seg(weight, L, R);
            }
        }
    }
}


void SegmentGraph::save_seg_pll_index(char *filename) {
    std::ofstream fout(filename, std::ios::binary);
    for (int i = 0; i < node_num_; i++) {
        unsigned size = hopIndex[i].size();
        fout.write((char *) &size, sizeof(unsigned));
        for (int j = 0; j < size; j++) {
            unsigned id = hopIndex[i][j].first;
            fout.write((char *) &id, sizeof(unsigned));
            unsigned list_size = hopIndex[i][j].second.size();
            fout.write((char *) &list_size, sizeof(unsigned));
            for (auto u: hopIndex[i][j].second.list) {
                fout.write((char *) &u.weight, sizeof(double));
                fout.write((char *) &u.L, sizeof(unsigned));
                fout.write((char *) &u.R, sizeof(unsigned));
            }
        }
    }
}


void SegmentGraph::save_segraph_index_add(char *filename) {
    std::ofstream fout(filename, std::ios::binary);
    unsigned check = 0, length = brute_tree_nodes.size();
    fout.write((char *) &length, sizeof(unsigned));
    std::cerr<<length<<endl;
    for (int i = 0; i < length; i++) {
        for (auto &node: brute_tree_nodes[i]) {;
            fout.write((char *) &node.height, sizeof(int));
            fout.write((char *) &node.father, sizeof(int));
            unsigned cnt = node.ch.size();
            check+= cnt;
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
        fout.write((char *) preserve_order[i].data(), sizeof(int) * node_num_);
    }
}


void SegmentGraph::load_segraph_index_add(char *filename) {
    std::ifstream fin(filename, std::ios::binary);
    unsigned length;
    fin.read((char *) &length, sizeof(unsigned));
    brute_tree_nodes.resize(length);
    preserve_order.resize(length);
    std::cerr<<length<<endl;
    for (unsigned cur = 0; cur < length; cur++) {
        std::cerr<<cur<<endl;
        unsigned cnt;
        unsigned check = 0;
        brute_tree_nodes[cur].resize(node_num_);
        for (auto &node: brute_tree_nodes[cur]) {
            check++;
            fin.read((char *) &node.height, sizeof(int));
            fin.read((char *) &node.father, sizeof(int));
            fin.read((char *) &cnt, sizeof(unsigned));
            node.ch.resize(cnt);
            node.ch.reserve(cnt);
            fin.read((char *) node.ch.data(), sizeof(unsigned) * cnt);
            unsigned edge_cnt;
            check+=cnt*4;
            fin.read((char *) &edge_cnt, sizeof(unsigned));
            for (int i = 0; i < edge_cnt; i++) {
                int id;
                fin.read((char *) &id, sizeof(int));
                unsigned list_cnt;
                fin.read((char *) &list_cnt, sizeof(unsigned));
                node.edges[id].list.reserve(list_cnt);
                for (int j = 0; j < list_cnt; j++) {
                    double weight;
                    unsigned L, R;
                    fin.read((char *) &weight, sizeof(double));
                    fin.read((char *) &L, sizeof(unsigned));
                    fin.read((char *) &R, sizeof(unsigned));
                    node.edges[id].emplace_back_seg(weight, L, R);
                }
            }
        }
        preserve_order[cur].resize(node_num_);
        fin.read((char *) preserve_order[cur].data(), sizeof(int) * node_num_);
    }
    segment_mapping(Left_Bound_, Right_Bound_);
}




