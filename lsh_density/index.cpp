#include "index.h"



void build_index(std::unordered_map<int, Traj>& traj_dict,
    int len, 
    int step, 
    double interval, 
    LSH_KDE_BKS& index)
{
    for(auto&& p:traj_dict){
        auto& t = p.second;
        auto& id = p.first;
        sliding_window(t, len, step, interval, [&](const TrajView& tv){
            int offset = std::distance(&t[0], tv.p);
            std::vector<double> flat_spatial_vec = tv.to_flat_vec_spatial();

            hash_combine(offset, id);
            index.insert(flat_spatial_vec, offset);
        });
    }
}



void build_index(std::unordered_map<int, Traj>& traj_dict,
    int len, 
    int step, 
    double interval, 
    LSH_KDE& index)
{
    for(auto&& p:traj_dict){
        auto& t = p.second;
        // auto& id = p.first;
        sliding_window(t, len, step, interval, [&](const TrajView& tv){
            std::vector<double> flat_spatial_vec = tv.to_flat_vec_spatial();
            index.insert(flat_spatial_vec);
        });
    }
}

