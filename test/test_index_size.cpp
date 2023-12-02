#include "Graph.h"
#include "SegmentGraph.h"
#include <getopt.h>
#include "utils.h"
#include "Segment.h"
#include <fstream>
#include <iostream>
#include <string>

using namespace std;

size_t getCurrentRSS() {
    std::ifstream statm("/proc/self/statm");
    size_t size, resident, share, text, lib, data, dt;
    statm >> size >> resident >> share >> text >> lib >> data >> dt;
    long page_size = sysconf(_SC_PAGESIZE);
    return resident * page_size;
}

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
    size_t size1 = getCurrentRSS();
    std::cerr << size1 << endl;
    if (method_type == 1) {
        graph.load_segraph_index(index_path);
    } else if (method_type == 2) {
        graph.load_segraph_index_add(index_path);
    } else if (method_type == 3) {
        graph.load_seg_pll_index(index_path);
    }
    size_t size2 = getCurrentRSS();
    std::cerr << size2 << endl;
    std::ofstream fout(save_log_path, ios::app);
    fout<<"current dataset:: "<<data_name<<endl;
    if (method_type == 1)
        fout << "NC index memory(MB) gap length: " << attr_length << " size->:: " << (size2 - size1) / (1024 * 1024) << endl;
    else if (method_type == 2)
        fout << "FUll index memory(MB) gap length: " << attr_length << " size->:: " << (size2 - size1) / (1024 * 1024) << endl;
    else if (method_type == 3)
        fout << "PLL memory(MB) gap length: " << attr_length << " size->:: " << (size2 - size1) / (1024 * 1024) << endl;

//    freopen(save_log_path, "a", stdout);
//    if (method_type == 1 || method_type == 2) {
//        printf("current data set %s\n", data_name);
//        printf("default attr length:: %d\n", attr_length);
//        if (method_type == 1)
//            printf("ave tree index search time(s) use::");
//        else
//            printf("ave full index search time(s) use::");
//    }
//    std::fclose(stdout);
    return 0;
}