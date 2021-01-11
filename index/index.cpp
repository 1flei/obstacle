#include "index.h"



void build_index(std::unordered_map<int, Traj>& traj_dict,
    int len, 
    int step, 
    double interval, 
    hnswlib::HierarchicalNSW<double, SubtIdx>& hnsw, 
    int NTHREADS)
{
    std::vector<std::thread> thread_pool;
    thread_pool.reserve(NTHREADS);

    for(int i=0;i<NTHREADS;i++){
        thread_pool.emplace_back([&, i, NTHREADS](){
            int cnt = 0;
            for(auto&& p:traj_dict){
                auto& t = p.second;
                auto& id = p.first;
                ++cnt;
                if(cnt%NTHREADS!=i){
                    continue;
                }
                // fmt::print("traj-{}, lens={}, #processed_traj={}\n", id, t.get_num_pnts(), cnt);
                sliding_window(t, len, step, interval, [&](const TrajView& tv){
                    int offset = std::distance(&t[0], tv.p);
                    hnsw.addPoint(tv.p, SubtIdx(id, offset));
                });
            }
        });
    }
    for(int i=0;i<NTHREADS;i++){
        thread_pool[i].join();
    }
}



SubtPriQueue searchKNN(hnswlib::HierarchicalNSW<double, SubtIdx>& hnsw, void* query, int k)
{
    SubtPriQueue myheap; //min-heap

    hnsw.forCandidates(query, k, [&](int , const SubtIdx& external_id, double dist) -> double{
        myheap.emplace_k(dist, external_id, std::max<int>(hnsw.ef_, k));
        if(myheap.size() >= k) {
            return myheap.maximum().first;
        }
        return std::numeric_limits<double>::max();
    });
    return myheap;
}