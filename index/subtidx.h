#pragma once

#include "../util.h"
#include <iostream>
#include <tuple>

struct SubtIdx
{
    int tid, offset;
    SubtIdx() : tid(0), offset(0) {}
    SubtIdx(int tid, int offset): tid(tid), offset(offset) {}

    bool operator<(const SubtIdx& idx) const {
        return tid < idx.tid || (tid==idx.tid && offset < idx.offset);
    }
    bool operator==(const SubtIdx& idx) const {
        return tid == idx.tid && offset == idx.offset;
    }
    bool operator!=(const SubtIdx& idx) const {
        return tid != idx.tid || offset != idx.offset;
    }
    template<class OS>
    friend OS& operator<<(OS& o, const SubtIdx& idx) {
        return o << "(" << idx.tid << "," << idx.offset << ")";
    }
};
//---------------------------------------------------------
//tuple hasher
// function has to live in the std namespace 
// so that it is picked up by argument-dependent name lookup (ADL).
namespace std{
    template<>
    struct hash<SubtIdx> 
    {
        size_t operator()(SubtIdx const& idx) const
        {                                              
            auto hasher = std::hash<std::tuple<int, int> >();                               
            return hasher(std::make_tuple(idx.tid, idx.offset)); 
        }
    };
}

//my wrapper to hnswlib
// class HNSW
// {
// public:
//     HNSW() {}
//     HNSW(const std::string& filename) {}
//     HNSW(const std::string& filename, const std::string& collectionName) {}


//     std::unique_ptr<hnswlib::HierarchicalNSW<double>> indexp;
// };

inline void sliding_window_trajdict(const std::unordered_map<int, Traj>& traj_dict, 
    int len, 
    int step, 
    double interval, 
    const tl::function_ref<void(const SubtIdx&, const TrajView&)>& f
)
{
    for(auto&& p:traj_dict){
        const Traj& t = p.second;
        int tid = p.first;
        sliding_window(t, len, step, interval, [&](const TrajView& tv){
            const Pnt3* p = tv.p;
            int offset = std::distance(&t[0], p);
            f(SubtIdx(tid, offset), tv);
        });
    }
}