#include "Graph.h"
#include "SegmentGraph.h"
#include <getopt.h>
#include "utils.h"
#include "Segment.h"

using namespace std;

int main(int argc, char *argv[]) {

    const struct option longopts[] = {
            // graph Path
            {"data_path",      required_argument, 0, 'g'},
            {"PLL_index_path", no_argument,       0, 'p'},
            {"CH_index_path",  no_argument,       0, 'c'},
            {"H2H_index_path", no_argument,       0, 'h'}
    };

    int ind;
    int iarg = 0;
    opterr = 1;    //getopt error message (off: 0)

    char PLL_index_path[256] = "";
    char CH_index_path[256] = "";
    char H2H_index_path[256] = "";
    char graph_path[256] = "";


    while (iarg != -1) {
        iarg = getopt_long(argc, argv, "g:p:c:h:", longopts, &ind);
        switch (iarg) {
            case 'g':
                if (optarg) {
                    strcpy(graph_path, optarg);
                }
                break;
            case 'p':
                if (optarg) {
                    strcpy(PLL_index_path, optarg);
                }
                break;
            case 'c':
                if (optarg) {
                    strcpy(CH_index_path, optarg);
                }
                break;
            case 'h':
                if (optarg) {
                    strcpy(H2H_index_path, optarg);
                }
                break;
        }
    }
    SegmentGraph graph;
    graph.load_weight_graph(graph_path);
    graph.add_random_range_attr(1, 10);
    if (isFileExists_ifstream(CH_index_path)) {
        graph.load_segraph_index(CH_index_path);
        std::cerr << "load graph finished" << endl;
    } else {
        graph.build_range_tree_decompose(2);
        graph.save_segraph_index(CH_index_path);
    }
    srand(0);
    std::vector<std::pair<unsigned, unsigned>> random_query;
    std::vector<double> result;
    for (int i = 0; i < 1000; i++) {
        unsigned u = rand() % graph.weight_graph_.size() + 1;
        unsigned v = rand() % graph.weight_graph_.size() + 1;
        random_query.emplace_back(u, v);
    }
    auto s = chrono::high_resolution_clock::now();
    for (int i = 0; i < 1000; i++) {
        unsigned u = random_query[i].first;
        unsigned v = random_query[i].second;
        Segment Seg = std::make_pair(1, 9);
        double ans1 = graph.dijkstra_range_path(u, v, graph.weight_graph_, Seg);
        if(ans1>=FLOAT_IVF) ans1 = -1;
        result.push_back(ans1);
    }
    auto e = chrono::high_resolution_clock::now();
    chrono::duration<double> diff = e - s;
    double time_slap = diff.count();
    std::cerr << "dijkstra time use:: " << time_slap << endl;

    s = chrono::high_resolution_clock::now();
    int cnt = 0;
    for (auto node_pair: random_query) {
        unsigned u = node_pair.first;
        unsigned v = node_pair.second;
        Segment Seg = std::make_pair(1, 9);
        double ans2 = graph.bi_dijkstra_range_path(u, v, graph.weight_graph_, Seg);
        if (ans2 >= FLOAT_IVF) ans2 = -1;
        if (!double_equal(result[cnt], ans2)) {
            std::cerr<<result[cnt]<<" "<<ans2<<endl;
            std::cerr << "bi dijkstra error" << endl;
            return 0;
        }
        cnt++;
    }
    e = chrono::high_resolution_clock::now();
    diff = e - s;
    time_slap = diff.count();
    std::cerr << "bi dijkstra time use:: " << time_slap << endl;

    s = chrono::high_resolution_clock::now();
    cnt = 0;
    for (auto node_pair: random_query) {
        unsigned u = node_pair.first;
        unsigned v = node_pair.second;
        double ans5 = graph.range_tree_decompose_query(u, v, 1, 9);
        if (ans5 >= FLOAT_IVF) ans5 = -1;
        if (!double_equal(result[cnt], ans5)) {
            std::cerr << "H2H error" << endl;
            return 0;
        }
        cnt++;
    }
    e = chrono::high_resolution_clock::now();
    diff = e - s;
    time_slap = diff.count();
    std::cerr << "Tree decomposed time use:: " << time_slap << endl;
    return 0;
}