#include "Graph.h"
#include "SegmentGraph.h"
#include <getopt.h>
#include "utils.h"
#include "Segment.h"
using namespace std;




int main(int argc, char *argv[]) {

    const struct option longopts[] = {
            // File Path
            {"data_path",              required_argument, 0, 'g'},
            {"Save_index_path",        no_argument,       0, 'i'},
            {"Save_log_path",          no_argument,       0, 's'},
            {"data_set_name",          no_argument,       0, 'n'},
            // Method Parameter
            {"attribute length",       required_argument, 0, 'l'},
            {"index construct method", required_argument, 0, 'm'},
    };

    int ind;
    int iarg = 0;
    opterr = 1;    //getopt error message (off: 0)
    unsigned attr_length, method_type;
    char index_path[256] = "";
    char save_log_path[256] = "";
    char graph_path[256] = "";
    char data_name[256] = "";


    while (iarg != -1) {
        iarg = getopt_long(argc, argv, "g:i:s:l:m:n:", longopts, &ind);
        switch (iarg) {
            case 'g':
                if (optarg) {
                    strcpy(graph_path, optarg);
                }
                break;
            case 'i':
                if (optarg) {
                    strcpy(index_path, optarg);
                }
                break;
            case 's':
                if (optarg) {
                    strcpy(save_log_path, optarg);
                }
                break;
            case 'n':
                if (optarg) {
                    strcpy(data_name, optarg);
                }
                break;
            case 'l':
                if (optarg) attr_length = atoi(optarg);
                break;
            case 'm':
                if (optarg) method_type = atoi(optarg);
                break;
        }
    }
    SegmentGraph graph;
    graph.load_weight_graph(graph_path);
    graph.add_random_range_attr(1, attr_length);
    if (method_type == 1) {
        graph.load_segraph_index(index_path);
    } else if (method_type == 2) {
        graph.load_segraph_index_add(index_path);
    } else if (method_type == 3) {
        graph.load_seg_pll_index(index_path);
    }

    std::vector<double> dijkstra_time, bi_dijkstra_time, tree_time, PLL_time;
    std::vector<pair<unsigned, unsigned >> query_dist;
    unsigned count_bound = 10000;
    dijkstra_time.resize(count_bound + 1);
    bi_dijkstra_time.resize(count_bound + 1);
    tree_time.resize(count_bound + 1);
    PLL_time.resize(count_bound + 1);
    query_dist.resize(count_bound + 1);
    std::cerr << "begin test::" << endl;

    srand(123);
    for (unsigned cur_count = 1; cur_count <= count_bound; cur_count++) {
        double ans1 = -1, time_slap;
        unsigned left, right, u, v;
        Segment Seg;
        if (cur_count % 100 == 0) std::cerr << cur_count << endl;

        unsigned length = (cur_count-1)/1000 + 1;
        graph.random_walk_query(u,v,left,right, length);
        auto s = chrono::high_resolution_clock::now();
        if (method_type == 1) {
            ans1 = graph.range_tree_decompose_query(u, v, left, right);
        } else if (method_type == 2) {
            ans1 = graph.brute_range_decompose_query(u, v, left, right);
        } else if (method_type == 3) {
            ans1 = graph.range_PLL_query(u, v, left, right);
        }
        if (ans1 >= FLOAT_IVF) ans1 = -1;
        auto e = chrono::high_resolution_clock::now();
        chrono::duration<double> diff = e - s;
        time_slap = diff.count();
        Seg = std::make_pair(left, right);
        tree_time[cur_count] = time_slap;

        double ans2 = -1;
        s = chrono::high_resolution_clock::now();
        ans2 = graph.bi_dijkstra_range_path(u, v, graph.weight_graph_, Seg);
        e = chrono::high_resolution_clock::now();
        diff = e - s;
        time_slap = diff.count();
        bi_dijkstra_time[cur_count] = time_slap;


        double ans3 = -1;
        s = chrono::high_resolution_clock::now();
        ans3 = graph.dijkstra_range_path(u, v, graph.weight_graph_, Seg);
        e = chrono::high_resolution_clock::now();
        diff = e - s;
        time_slap = diff.count();
        dijkstra_time[cur_count] = time_slap;
        if (!double_equal(ans1, ans2))
            std::cerr << "ERROR:: " << u << " " << v << " " << ans1 << " " << ans2 << endl;
        if (!double_equal(ans2, ans3))
            std::cerr << "ERROR:: " << u << " " << v << " " << ans1 << " " << ans2 << endl;
        query_dist[cur_count].first = length;
        query_dist[cur_count].second = cur_count;
    }

    double time_dij = 0.0, time_bi_dij = 0.0, time_tree = 0.0;
    for (unsigned cur = 1; cur <= count_bound; cur++) {
        time_dij += dijkstra_time[cur];
        time_bi_dij += bi_dijkstra_time[cur];
        time_tree += tree_time[cur];
    }
    time_tree /= count_bound;
    time_bi_dij /= count_bound;
    time_dij /= count_bound;

    std::vector<double> split_time_tree(11), split_time_bi_dij(11), split_time_dij(11);

    unsigned cur = 1;
    sort(query_dist.begin(), query_dist.end());
    int length_tag[100] = {0};
    for (int i = 1; i <= 10; i++) {
        for (int j = 1; j <= count_bound / 10; j++) {
            unsigned idx = query_dist[cur].second;
            unsigned length = query_dist[cur].first;
            length_tag[query_dist[cur].first]++;
            split_time_tree[length] += tree_time[idx];
            split_time_bi_dij[length] += bi_dijkstra_time[idx];
            split_time_dij[length] += dijkstra_time[idx];
            cur++;
        }
        std::cerr << query_dist[cur].first << endl;
    }
    for (int i = 1; i <= 10; i++) {
        if (length_tag[i]) {
            split_time_tree[i] /= length_tag[i];
            split_time_bi_dij[i] /= length_tag[i];
            split_time_dij[i] /= length_tag[i];
        }
    }
    freopen(save_log_path, "a", stdout);
    if (method_type == 1 || method_type == 2) {
        printf("current data set %s\n", data_name);
        printf("default attr length:: %d\n", attr_length);
        if (method_type == 1)
            printf("ave tree index search time(s) use:: %.6f\n", time_tree);
        else
            printf("ave full index search time(s) use:: %.6f\n", time_tree);

        printf("tree index search time(s): sort by distance split:10\n");
        for (int i = 1; i <= 10; i++) printf("%.6f ", split_time_tree[i]);
        printf("\n");
        printf("ave bi-dij time(s) use:: %.6f\n", time_bi_dij);
        printf("bi-dij search time(s): sort by distance split:10\n");
        for (int i = 1; i <= 10; i++) printf("%.6f ", split_time_bi_dij[i]);
        printf("\n");
        printf("ave dij time(s) use:: %.6f\n", time_dij);
        printf("dij search time(s): sort by distance split:10\n");
        for (int i = 1; i <= 10; i++) printf("%.6f ", split_time_dij[i]);
        printf("\n");
        printf("*************************\n");
    } else {
        printf("current data set %s\n", data_name);
        printf("default attr length:: %d\n", attr_length);
        printf("ave PLL index search time(s) use:: %.6f\n", time_tree);

        printf("PLL index search time(s): sort by distance split:10\n");
        for (int i = 1; i <= 10; i++) printf("%.6f ", split_time_tree[i]);
        printf("\n");
        printf("ave bi-dij time(s) use:: %.6f\n", time_bi_dij);
        printf("bi-dij search time(s): sort by distance split:10\n");
        for (int i = 1; i <= 10; i++) printf("%.6f ", split_time_bi_dij[i]);
        printf("\n");
        printf("ave dij time(s) use:: %.6f\n", time_dij);
        printf("dij search time(s): sort by distance split:10\n");
        for (int i = 1; i <= 10; i++) printf("%.6f ", split_time_dij[i]);
        printf("\n");
        printf("*************************\n");
    }
    std::fclose(stdout);
    return 0;
}