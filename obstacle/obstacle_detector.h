#pragma once

#include <unordered_map>
#include "../util.h"
#include "lsh_density_patch.h"
#include "../index/index.h"
#include <thread>
#include "json.hpp"
#include "../my_timer.h"
#include <stack>

class ObstacleDetector
{
public:
    using json = nlohmann::json;
public:
    ObstacleDetector(
        int lsh_density_nrepeat, 
        double sliding_window_length, 
        double sliding_window_step, 
        double sliding_window_interval, 
        double sigma, 
        double min_support, 
        double significance, 
        int hnsw_max_objs, 
        int hnsw_m, 
        int hnsw_ef_construction,
        int hnsw_build_n_threads = 16, 
        int num_neighbors = 32
    );
    ~ObstacleDetector();

    //build index for density query given a reference dataset
    //will build lsh_density for reference set, reference set succeed sub-trajectory, query set, and query set succeed sub-trajectory
    //then build hnsw 

    void build_lsh_density(LSH_KDE& lsh_kde, LSH_KDE& lsh_kde_succ, const std::unordered_map<int, Traj>& traj_dict, int& n_cnt);
    void build_lsh_density_ref(const std::unordered_map<int, Traj>& traj_dict) {
        build_lsh_density(lsh_kde_ref, lsh_kde_succ_ref, traj_dict, lsh_kde_n_ref);
    }
    void build_lsh_density_q(const std::unordered_map<int, Traj>& traj_dict) {
        build_lsh_density(lsh_kde_q, lsh_kde_succ_q, traj_dict, lsh_kde_n_q);
    }

    void build_hnsw(hnswlib::HierarchicalNSW<double, SubtIdx>& hnsw, const std::unordered_map<int, Traj>& traj_dict);

    void build_hnsw_ref(const std::unordered_map<int, Traj>& traj_dict) {
        build_hnsw(hnsw_ref, traj_dict);
    }
    void build_hnsw_q(const std::unordered_map<int, Traj>& traj_dict) {
        build_hnsw(hnsw_q, traj_dict);
    }

    void build_index_ref(const std::unordered_map<int, Traj>& traj_dict) {
        ref_dict_ptr = &traj_dict;
        fmt::print("building lsh_kde\n");
        build_lsh_density_ref(traj_dict);
        fmt::print("building hnsw\n");
        build_hnsw_ref(traj_dict);
    }
    void build_index_q(const std::unordered_map<int, Traj>& traj_dict) {
        queries_dict_ptr = &traj_dict;
        build_lsh_density_q(traj_dict);
        // build_hnsw_q(traj_dict);
    }

    void save_index(const std::string& filename);
    void save_index_q(const std::string& filename);
    void load_index(const std::string& filename, const std::unordered_map<int, Traj>& traj_dict);
    void load_index_q(const std::string& filename, const std::unordered_map<int, Traj>& traj_dict);


    //query obstacle given the index and a trajectory set of interest
    void detect_obstacle(const std::unordered_map<int, Traj>& query_traj_dict);

    inline void for_consecutive_subtraj_spatial_vec(
        const Traj& traj, 
        const tl::function_ref<void(const std::vector<double>&, const std::vector<double>&)>& f
    ) 
    {
        std::vector<double> cur_traj_pnts;
        std::vector<double> prev_traj_pnts;
        int prev_offset = -100;
        enum_sliding_window(traj, sliding_window_length, sliding_window_step, sliding_window_interval, [&](int offset, const TrajView& tv){
            //will copy 
            cur_traj_pnts = tv.to_flat_vec_spatial();
            // fmt::print("offset={}, prev_offset={}\n", offset, prev_offset);
            if(offset <= prev_offset+2){
                f(prev_traj_pnts, cur_traj_pnts);
            } else if(prev_offset>=0){
                f(prev_traj_pnts, prev_traj_pnts);
            }

            prev_offset = offset;
            std::swap(prev_traj_pnts, cur_traj_pnts);
        });
        //last subt
        if(prev_offset>=0){
            f(prev_traj_pnts, prev_traj_pnts);
        }
    }

    
    inline void sliding_window_trajdict(const std::unordered_map<int, Traj>& traj_dict, 
        tl::function_ref<void(const SubtIdx&, const TrajView&)> f)
    {
        for(auto&& p:traj_dict){
            const Traj& t = p.second;
            int tid = p.first;
            sliding_window(t, sliding_window_length, sliding_window_step, sliding_window_interval, [&](const TrajView& tv){
                const Pnt3* p = tv.p;
                int offset = std::distance(&t[0], p);
                f(SubtIdx(tid, offset), tv);
            });
        }
    }

    double get_sigma() {
        return sigma;
    }

    int get_sliding_window_length() {
        return sliding_window_length;
    }

    double get_density_ref(const std::vector<double>& q) {
        return lsh_kde_ref.get_estimated_density(q);
    }

    double get_density_succ_ref(const std::vector<double>& q, const std::vector<double>& q_succ) {
        return get_estimated_succ_density(lsh_kde_succ_ref, q, q_succ);
    }

    double get_density_q(const std::vector<double>& q) {
        return lsh_kde_q.get_estimated_density(q);
    }

    double get_density_succ_q(const std::vector<double>& q, const std::vector<double>& q_succ) {
        return get_estimated_succ_density(lsh_kde_succ_q, q, q_succ);
    }

    void get_obstacle_one_subt(
        const SubtIdx& qidx
        );

        
    nlohmann::json output_obstacle_json(const std::unordered_map<int, Traj>& query_dict);
    nlohmann::json output_obstacle_json_with_index(const std::unordered_map<int, Traj>& query_dict, const std::string& filename);

    std::tuple<TrajView, TrajView> get_view_from_subtidx(const std::unordered_map<int, Traj>& traj_dict, const SubtIdx& idx);

    //get f(q), f_succ(q) for ref and q given two sub-trajectory views 
    std::tuple<double, double, double, double> get_densities(
        const TrajView& tv_nr, 
        const TrajView& tv_nr_succ
    );

    bool obstacle_conditaion(double dq, double dr, double dq_comover, double dr_comover)
    {
        // double ratio = dr_comover*lsh_kde_n_ref/dr/dr / (dq_comover*lsh_kde_n_q/dq/dq);
        double ratio = (dr_comover/dr) / (dq_comover/dq);
        return ratio > significance && dq > min_support && dr_comover > min_support;
    }

    
    nlohmann::json get_single_obstacle_json(
        double dq, 
        double dr, 
        double dq_comover, 
        double dr_comover, 
        const TrajView& tv_nr_succ
    )
    {
        nlohmann::json j;
        Pnt3 lastp = tv_nr_succ[tv_nr_succ.get_num_pnts()-1];
        std::vector<double> tmp{lastp.t, lastp.x, lastp.y};
        // fmt::print("{}, {}, {}\n", lastp.t, lastp.x, lastp.y);
        j["obs_pnts"] = tmp;
        j["t_succ"] = tv_nr_succ.to_flat_vec();

        
        std::vector<double> densities{dq, dr, dq_comover, dr_comover};
        j["densities"] = densities;
        return j;
    }

    int get_num_ref_pnts() {
        return lsh_kde_n_ref;
    }
    int get_num_q_pnts() {
        return lsh_kde_n_q;
    }

    //let a be root
    void merge_obstacle(const SubtIdx& a, const SubtIdx& b) {
        auto pa = obst_subtidx_belongs_to.get(a);
        auto pb = obst_subtidx_belongs_to.get(b);
        if(pa==pb){
            return;
        }

        obst_subtidx_belongs_to.merge(a, b);
        auto ita = obst_json_map.find(pa);
        auto itb = obst_json_map.find(pb);

        if(itb!=obst_json_map.end()) {
            if(ita!=obst_json_map.end()){ 
                //both non-empty
                // fmt::print("a.size()=={};  b.size()=={}\n", ita->second.size(), itb->second.size());
                if(ita->second.size() < itb->second.size()){ 
                    std::swap(ita->second, itb->second);
                }
                for(auto&& j:itb->second){ 
                    ita->second.push_back(std::move(j));
                }
                obst_json_map.erase(itb);
            } else{
                obst_json_map[a] = std::move(itb->second);
                obst_json_map.erase(itb);
            }
        }
    }

    LSH_KDE& get_lsh_kde_ref() {
        return lsh_kde_ref;
    }

private:
    constexpr static double lsh_density_normal_r = 5.12;
    const static int lsh_density_normal_k = 8;
    const static int lsh_density_normal_l = 3;

    size_t dim;
    double sliding_window_length, sliding_window_step, sliding_window_interval;
    double sigma;
    double min_support;
    double significance;
    int hnsw_build_n_threads;
    int num_neighbors;
    double dist2_threshold;
    int hnsw_max_objs;

    int lsh_kde_n_ref;
    int lsh_kde_n_q;

    Pnt3EuSpace hnsw_space;

    const std::unordered_map<int, Traj>* ref_dict_ptr;
    const std::unordered_map<int, Traj>* queries_dict_ptr; 

    LSH_KDE lsh_kde_ref;
    LSH_KDE lsh_kde_succ_ref;

    LSH_KDE lsh_kde_q;
    LSH_KDE lsh_kde_succ_q;
    
    hnswlib::HierarchicalNSW<double, SubtIdx> hnsw_ref; 
    hnswlib::HierarchicalNSW<double, SubtIdx> hnsw_q; 

    std::unordered_set<SubtIdx> q_checked;
    std::unordered_map<SubtIdx, bool> ref_is_candidates; 

    std::unordered_map<SubtIdx, json> obst_json_map;
    HashDisjointSet<SubtIdx> obst_subtidx_belongs_to;
};