//
// Created by Mingyu on 23-10-20.
//

#include "SegmentGraph.h"
#include "Graph.h"
#include "utils.h"

using namespace std;

/*
 * use degree decent order
 */

void SegmentGraph::build_range_PLL_index(unsigned left_range, unsigned right_range) {
    vector<pair<unsigned, unsigned >> deg2corevid;
    unordered_map<unsigned , unsigned > corevid2level;
    vector<unsigned > level2corevid;
    for (int i = 0; i < node_num_; i++) {
        deg2corevid.emplace_back(weight_graph_[i].size(), i);
    }
    sort(deg2corevid.begin(), deg2corevid.end(), std::greater<pair<unsigned, unsigned>>());
    for (auto &p: deg2corevid) level2corevid.push_back(p.second);

    for (unsigned level = 0; level < deg2corevid.size(); ++level)corevid2level[deg2corevid[level].second] = level;


    for (int i = 0; i < node_num_; i++) hopIndex[i];
    unsigned check_tag = 0;
    for (auto &covid: level2corevid) {
        check_tag++;
        unsigned level = corevid2level[covid];
        std::cerr<<check_tag<<endl;
        if (check_tag % (node_num_ / 1000 + 1) == 0) std::cerr << "PLL current:: " << check_tag <<endl;
        auto &sidx = hopIndex[covid];
        vector<priority_queue<SegQueNode>> PrioQs(right_range - left_range + 1 + 1);
        unsigned lsize = 0;
        PrioQs[lsize].emplace(covid, right_range, left_range, 0);
        while (!PrioQs[lsize].empty()) {
            SegQueNode v = PrioQs[lsize].top();
            PrioQs[lsize].pop();
            auto &didx = hopIndex[v.id];
            if (!dominate_check(sidx, didx, v)) {//not dominated, insert
                if (didx.empty() || didx.back().first != level) {
                    didx.emplace_back(level, SegList());
                }

                didx.back().second.emplace_back_seg(v.weight, v.L, v.R);

                for (unsigned j = 0; j < weight_graph_[v.id].size(); j++) {
                    auto e = weight_graph_[v.id][j];
                    unsigned attr = attr_graph_[v.id][j];
                    if (corevid2level[e.first] <= level) continue;
                    double weight = v.weight + e.second;
                    unsigned l = std::min(attr, v.L);
                    unsigned r = std::max(attr, v.R);
                    PrioQs[r - l + 1].emplace(e.first, l, r, weight);
                }
            }
            while (lsize < PrioQs.size() && PrioQs[lsize].empty()) lsize++;
            if (lsize == PrioQs.size())break;
        }
    }
    std::cerr<<"PLL range finished"<<endl;
}