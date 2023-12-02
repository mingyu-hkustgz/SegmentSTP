//
// Created by BLD on 23-10-19.
//
#include "utils.h"

#ifndef SEGMENTSTP_SEGMENT_H
#define SEGMENTSTP_SEGMENT_H
using namespace std;

struct Sedge {

    Sedge(double w, unsigned l, unsigned r) {
        weight = w;
        L = l;
        R = r;
    }


    double weight;

    unsigned L, R;

    bool operator<(const Sedge &a) const {
        return this->weight < a.weight;
    }

    inline Sedge operator+(const Sedge &a) const {
        Sedge res{this->weight + a.weight, std::min(a.L, this->L), std::max(a.R, this->R)};
        return res;
    }

    bool check_in_segment(unsigned left_range, unsigned right_range) const;

    bool check_be_cover(Sedge &a) const;

};

struct SegList {

    inline SegList &emplace_back_seg(double weight, unsigned L, unsigned R) {
        list.emplace_back(weight, L, R);
        return *this;
    }

    inline bool operator==(const SegList &b) {
        if (list.size() != b.size()) return false;
        for (unsigned i = 0; i < list.size(); i++) {
            if (list[i].L != b.list[i].L) return false;
            if (list[i].R != b.list[i].R) return false;
            if (!double_equal(list[i].weight, b.list[i].weight)) return false;
        }
        return true;
    }

    inline void combine(SegList &&a) {
        auto &la = list;
        auto &ra = a.list;
        std::vector<bool> lremoved(la.size(), false), rremoved(ra.size(), false);
        for (uint l = 0; l < la.size(); ++l) {
            for (uint r = 0; r < ra.size(); ++r) {
                if (rremoved[r])continue;
                if (la[l].weight < ra[r].weight) {
                    rremoved[r] = ra[r].check_be_cover(la[l]);
                } else if (la[l].weight > ra[r].weight) {
                    lremoved[l] = la[l].check_be_cover(ra[r]);
                } else {
                    rremoved[r] = ra[r].check_be_cover(la[l]);
                    if (!rremoved[r])
                        lremoved[l] = la[l].check_be_cover(ra[r]);
                }

                if (lremoved[l])break;
            }
        }
        std::vector<Sedge> tmp;
        uint l = 0, r = 0;
        while (l < la.size() || r < ra.size()) {
            if (l == la.size()) {
                if (!rremoved[r])
                    tmp.emplace_back(ra[r]);
                r++;
                continue;
            }
            if (r == ra.size()) {
                if (!lremoved[l])
                    tmp.emplace_back(la[l]);
                l++;
                continue;
            }
            if (la[l].weight <= ra[r].weight) {
                if (!lremoved[l])tmp.emplace_back(la[l]);
                l++;
            } else {
                if (!rremoved[r])tmp.emplace_back(ra[r]);
                r++;
            }
        }
        la = std::move(tmp);
    }

    inline SegList operator+(const SegList &a) const {
        SegList res;
        auto &o1 = this->list, &o2 = a.list;
        auto &r = res.list;
        for (auto &p1: o1) {
            for (auto &p2: o2) {
                res.emplace_back_seg(p1.weight + p2.weight, min(p1.L, p2.L), max(p1.R, p2.R));
            }
        }
        std::sort(r.begin(), r.end());
        res.redundant_reduction();
        return res;
    }

    void redundant_reduction() {
        std::vector<bool> removed(list.size(), false);
        for (unsigned i = 0; i < list.size(); i++) {
            if (removed[i]) continue;
            for (unsigned j = i + 1; j < list.size(); j++) {
                if (removed[j]) continue;
                if (list[j].check_be_cover(list[i])) {
                    removed[j] = true;
                }
            }
        }
        unsigned p = 0;
        for (unsigned i = 0; i < list.size(); i++) {
            if (!removed[i]) {
                if (p < i)
                    list[p] = list[i];
                p++;
            }
        }
        list.erase(list.begin() + p, list.end());
    }

    double dist(unsigned left_range, unsigned right_range) {
        for (auto &p: list) {
            if (p.check_in_segment(left_range, right_range)) {
                return p.weight;
            }
        }
        return FLOAT_IVF;
    }


    inline unsigned size() const {
        return list.size();
    }

    std::vector<Sedge> list;
};

typedef std::vector<std::vector<unsigned> > AttrGraph;
typedef std::unordered_map<unsigned, SegList> SegAttr;
typedef std::vector<SegAttr> SegHop;

struct SegTreeNode {
    int height, father;
    std::vector<unsigned> ch;
    SegAttr edges;
};

struct SegQueNode {
    SegQueNode(unsigned i, unsigned l, unsigned r, double w) {
        id = i;
        L = l;
        R = r;
        weight = w;
    }

    unsigned id, L, R;
    double weight;

    inline bool operator<(const SegQueNode &u) const {
        return weight > u.weight;
    }
};

inline static bool dominate_check(vector<pair<unsigned int, SegList>> &s, vector<pair<unsigned int, SegList>> &d, const SegQueNode &node) {
    auto sit = s.begin(), dit = d.begin();
    while (sit != s.end() && dit != d.end()) {
        if (sit->first < dit->first)++sit;
        else if (dit->first < sit->first)++dit;
        else {
            for (auto &sp: sit->second.list) {
                if (sp.weight < node.weight) {
                    for (auto &dp: dit->second.list) {
                        if ((sp.weight + dp.weight) <= node.weight &&
                            node.L <= std::min(sp.L, dp.L) && node.R >= std::max(sp.R, dp.R))
                            return true;
                    }
                }
            }
            ++sit, ++dit;
        }
    }
    return false;
}

#endif //SEGMENTSTP_SEGMENT_H
