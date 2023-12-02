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
    char save_path[256] = "";
    unsigned part=1;

    while (iarg != -1) {
        iarg = getopt_long(argc, argv, "g:p:c:h:s:", longopts, &ind);
        switch (iarg) {
            case 'g':
                if (optarg) {
                    strcpy(graph_path, optarg);
                }
                break;
            case 'p':
                if (optarg) part = atoi(optarg);
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
            case 's':
                if (optarg) {
                    strcpy(save_path, optarg);
                }
                break;
        }
    }
    SegmentGraph graph;
    graph.load_weight_graph(graph_path);
    std::ofstream fout(save_path);
    unsigned range = graph.weight_graph_.size()*part/10;
    WeightGraph range_graph;
    range_graph.resize(range+1);
    unsigned edge_count =0;
    for(int i=0;i<range;i++){
        for(auto & j : graph.weight_graph_[i]){
            if(j.first<range){
                range_graph[i].push_back(j);
                edge_count++;
            }
        }
    }
    fout<<range<<" "<<edge_count<<endl;
    for(int i=0;i<range;i++){
        for(auto &j:range_graph[i]){
            fout<<i<<" "<<j.first<<" "<<j.second<<endl;
        }
    }

    return 0;
}