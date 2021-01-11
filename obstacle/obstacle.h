#pragma once

#include <vector>
#include <unordered_map>
#include "../index/index.h"
#include "json.hpp"
#include "../my_timer.h"

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <fstream>

inline void to_json(nlohmann::json& j, const Pnt3& p) {
    std::vector<double> tmp{p.t, p.x, p.y};
    j = tmp;
}
inline void from_json(const nlohmann::json& j, Pnt3& p) {
    p.t = j.at(0);
    p.x = j.at(1);
    p.y = j.at(2);
}

namespace boost {
namespace serialization {

template<class Archive>
void serialize(Archive & ar, SubtIdx& idx, const unsigned int )
{
    ar & idx.tid;
    ar & idx.offset;
}

}
}

// const double SIGMA = 0.1;
class ObsDetector
{
public:
    // const double SIGMA = 0.02;

    ObsDetector(
        std::unordered_map<int, Traj>& ref_dict, 
        std::unordered_map<int, Traj>& queries_dict, 
        hnswlib::HierarchicalNSW<double, SubtIdx>& hnsw, 
        hnswlib::HierarchicalNSW<double, SubtIdx>& hnsw_q, 
        hnswlib::SpaceInterface<double>& dtwspace, 
        int numNeighbors, 
        int verbosity = 0, 
        int swlen = 6, 
        int swstep = 1, 
        double swinterval = 600., 
        double sigma = 0.3, 
        double significance = 1.68, 
        double min_support = 2., 
        const std::string& ref_neighbor_cache_name = ""
    ) 
    : ref_dict(ref_dict), queries_dict(queries_dict), hnsw(hnsw), hnsw_q(hnsw_q), dtwspace(dtwspace), numNeighbors(numNeighbors), 
        verbosity(verbosity), swlen(swlen), swstep(swstep), swinterval(swinterval), sigma(sigma), dist2_threshold(sigma*sigma*9), 
        significance(significance), min_support(min_support), ref_neighbor_cache_name(ref_neighbor_cache_name)
    {
        if(ref_neighbor_cache_name!=""){
            std::ifstream ifs(ref_neighbor_cache_name);
            if(ifs.is_open()){
                boost::archive::binary_iarchive ia(ifs);

                ia & ref_neighbor_cache;
            }
        }
    }

    ~ObsDetector() {
        if(ref_neighbor_cache_name!=""){
            std::ofstream ofs(ref_neighbor_cache_name);
            if(ofs.is_open()){
                boost::archive::binary_oarchive oa(ofs);

                oa & ref_neighbor_cache;
            }
        }
    }


    std::vector<std::vector<Pnt3> > get_obstacles();
    
    //compute the obstacle point associated with qidx and store them in obstacle_pnts
    void get_obstacle_one_subt(
        const SubtIdx& qidx, 
        const TrajView& tv,  
        nlohmann::json& obstacle_pnts);

    nlohmann::json output_obstacle_json();

    void get_obstacle_one_subt_naive(
        const SubtIdx& qidx, 
        const TrajView& tv,  
        nlohmann::json& obstacle_pnts);

    nlohmann::json output_obstacle_json_naive();
    nlohmann::json output_subtraj_density_json_bbox(const BBox& bb);

    nlohmann::json get_single_obstacle_json(
        const TrajView& tv_nr, 
        const TrajView& tv_nr_succ, 
        double res_p_dist, 
        const SubtIdx& nr_idx, 
        const SubtIdx& qidx, 
        double dq, 
        double dr, 
        double dq_comover, 
        double dr_comover, 
        const SubtPriQueue& nqr, 
        const SubtPriQueue& nrr);

        

    inline void sliding_window_trajdict(std::unordered_map<int, Traj>& traj_dict, 
        tl::function_ref<void(const SubtIdx&, const TrajView&)> f)
    {
        for(auto&& p:traj_dict){
            Traj& t = p.second;
            int tid = p.first;
            sliding_window(t, swlen, swstep, swinterval, [&](const TrajView& tv){
                Pnt3* p = tv.p;
                int offset = std::distance(&t[0], p);
                f(SubtIdx(tid, offset), tv);
            });
        }
    }

    
    double density(const SubtPriQueue& res_que);
    double density_prev(std::unordered_map<int, Traj>& traj_dict, const TrajView& tv_prev, 
        const SubtPriQueue& res_que, hnswlib::SpaceInterface<double>& distf);
    double density_succ(std::unordered_map<int, Traj>& traj_dict, const TrajView& tv_succ, 
        const SubtPriQueue& res_que, hnswlib::SpaceInterface<double>& distf);
        
    double density_comover(std::unordered_map<int, Traj>& traj_dict, const TrajView& tv_succ, 
        const SubtPriQueue& res_que, hnswlib::SpaceInterface<double>& distf);


    double get_score(double dq, double dr, double dq_succ, double dr_succ);

    double get_score_fisher(double dq, double dr, double dq_comover, double dr_comover);

    double get_score_z(double dq, double dr, double dq_comover, double dr_comover);

    std::tuple<TrajView, TrajView> get_view_from_subtidx(std::unordered_map<int, Traj>& traj_dict, const SubtIdx& idx);

    std::tuple<double, double> get_densities(
        std::unordered_map<int, Traj>& traj_dict, 
        const TrajView& tv_succ, 
        const SubtPriQueue& res_que, 
        hnswlib::SpaceInterface<double>& distf,
        double dist_sigma
    );

    std::tuple<double, double, double, double> get_densities(
        const TrajView& tv_nr, 
        const TrajView& tv_nr_succ, 
        SubtPriQueue& nqr, 
        SubtPriQueue& nrr, 
        std::unordered_map<int, Traj>& ref_dict, 
        std::unordered_map<int, Traj>& queries_dict, 
        hnswlib::SpaceInterface<double>& dtwspace
    );

    void filter_res_que(
        std::unordered_map<int, Traj>& traj_dict, 
        SubtPriQueue& res_que
    );

    bool obstacle_conditaion(double dq, double dr, double dq_comover, double dr_comover)
    {
        double z = get_score_z(dq, dr, dq_comover, dr_comover);
        return z > significance && dq > min_support && dr > min_support;
    }


private:
    std::unordered_map<int, Traj>& ref_dict;
    std::unordered_map<int, Traj>& queries_dict; 
    hnswlib::HierarchicalNSW<double, SubtIdx>& hnsw; 
    hnswlib::HierarchicalNSW<double, SubtIdx>& hnsw_q; 
    hnswlib::SpaceInterface<double>& dtwspace; 
    int numNeighbors;
    int verbosity;
    
    std::unordered_map<SubtIdx, bool> ref_is_candidates; 
    std::unordered_set<SubtIdx> q_checked;

    std::unordered_map<SubtIdx, std::vector<std::pair<double, SubtIdx> > > ref_neighbor_cache;
    std::string ref_neighbor_cache_name;

    int swlen;
    int swstep;
    double swinterval;
    double sigma;
    double dist2_threshold;
    double significance; 
    double min_support;
};

