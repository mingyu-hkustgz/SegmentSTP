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
    unsigned attr_length, method_type, thread_num = 40;
    char index_path[256] = "";
    char save_log_path[256] = "";
    char graph_path[256] = "";
    char data_name[256] = "";


    while (iarg != -1) {
        iarg = getopt_long(argc, argv, "g:i:s:l:m:n:p:", longopts, &ind);
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
            case 'p':
                if (optarg) thread_num = atoi(optarg);
                break;
        }
    }
    SegmentGraph graph;
    graph.load_weight_graph(graph_path);
    graph.add_random_range_attr(1, attr_length);

    auto s = chrono::high_resolution_clock::now();
    if (method_type == 1) {
        graph.build_range_PLL_index(1, attr_length);
        graph.save_seg_pll_index(index_path);
    } else if (method_type == 2) {
        graph.build_range_tree_decompose(2);
        graph.save_segraph_index(index_path);
    } else if (method_type == 3) {
        graph.num_threads = thread_num;
        graph.build_range_tree_decompose(3);
    }
    else if (method_type == 4) {
        graph.num_threads = thread_num;
        graph.build_range_tree_decompose(4);
    }
    else if (method_type == 5) {
        graph.build_range_tree_decompose(5);
        graph.save_segraph_index_add(index_path);
    }else if (method_type == 6) {
        graph.build_brute_range_tree_and_save(index_path);
    }

    auto e = chrono::high_resolution_clock::now();
    freopen(save_log_path, "a", stdout);
    if (method_type >= 2) {
        std::vector<double> res;
        res = graph.tree_statistics_analysis();
        chrono::duration<double> diff = e - s;
        double time_use = diff.count();
        printf("current data set %s\n", data_name);
        printf("time use:: %.3f\n", time_use);
        printf("max_tree_width :: %.0f\n", res[0]);
        printf("sum_tree_width :: %.0f\n", res[1]);
        printf("max_tree_height :: %.0f\n", res[2]);
        printf("sum_tree_height :: %.0f\n", res[3]);
        printf("max_seg_size :: %.0f\n", res[4]);
        printf("ave_seg_size :: %.5f\n", res[5]);
    }else if(method_type==1){
        chrono::duration<double> diff = e - s;
        double time_use = diff.count();
        printf("current data set %s\n", data_name);
        printf("time use:: %.3f\n", time_use);
    }
    std::fclose(stdout);
    return 0;
}