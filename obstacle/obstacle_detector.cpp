#include "obstacle_detector.h"


ObstacleDetector::ObstacleDetector(
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
    int hnsw_build_n_threads, 
    int num_neighbors
)
    :dim(sliding_window_length*2), sliding_window_length(sliding_window_length), sliding_window_step(sliding_window_step), sliding_window_interval(sliding_window_interval),
        sigma(sigma), min_support(min_support), significance(significance), 
         hnsw_build_n_threads(hnsw_build_n_threads), dist2_threshold(0.25*sigma*sigma), lsh_kde_n_ref(0), lsh_kde_n_q(0), num_neighbors(num_neighbors), 
         hnsw_max_objs(hnsw_max_objs), 
        lsh_kde_ref(sliding_window_length*2, lsh_density_normal_k, lsh_density_normal_l, lsh_density_normal_r*sigma, lsh_density_nrepeat), 
        lsh_kde_succ_ref(sliding_window_length*2, lsh_density_normal_k, lsh_density_normal_l, lsh_density_normal_r*sigma, lsh_density_nrepeat), 
        lsh_kde_q(sliding_window_length*2, lsh_density_normal_k, lsh_density_normal_l, lsh_density_normal_r*sigma, lsh_density_nrepeat), 
        lsh_kde_succ_q(sliding_window_length*2, lsh_density_normal_k, lsh_density_normal_l, lsh_density_normal_r*sigma, lsh_density_nrepeat), 
        hnsw_space(sliding_window_length), 
        hnsw_ref(&hnsw_space, hnsw_max_objs, hnsw_m, hnsw_ef_construction),
        hnsw_q(&hnsw_space, hnsw_max_objs, hnsw_m, hnsw_ef_construction)
{
}

//build lsh_kde index for density query given a reference dataset
void ObstacleDetector::build_lsh_density(LSH_KDE& lsh_kde, LSH_KDE& lsh_kde_succ,const std::unordered_map<int, Traj>& traj_dict, 
    int& n_cnt)
{
    n_cnt = 0;
    for(auto&& p:traj_dict){
        int tid = p.first;
        const auto& traj = p.second;

        for_consecutive_subtraj_spatial_vec(traj, [&](const std::vector<double>& x, const std::vector<double>& x_succ) {
            lsh_kde.insert(x);
            insert_succ_density(lsh_kde_succ, x, x_succ);
            ++n_cnt;
        });
    }

    fmt::print("cnt={}\n", n_cnt);

    //now lsh_kde_ref is done, loop the trajview again
    
}

std::tuple<TrajView, TrajView> ObstacleDetector::get_view_from_subtidx(const std::unordered_map<int, Traj>& traj_dict, const SubtIdx& idx)
{
    auto it = traj_dict.find(idx.tid);
    assert(it != traj_dict.end());
    auto& traj_i = it->second;
    TrajView tv_nr(traj_i, idx.offset, sliding_window_length);
    TrajView tv_nr_succ = try_get_trajview_clamp01(traj_i, idx.offset+sliding_window_step, sliding_window_length, sliding_window_interval).value_or(tv_nr);
    return std::make_tuple(tv_nr, tv_nr_succ);
}


void ObstacleDetector::build_hnsw(
    hnswlib::HierarchicalNSW<double, SubtIdx>& hnsw, 
    const std::unordered_map<int, Traj>& traj_dict)
{
    std::vector<std::thread> thread_pool;
    thread_pool.reserve(hnsw_build_n_threads);

    int tot_cnt = 0;
    for(int i=0;i<hnsw_build_n_threads;i++){
        thread_pool.emplace_back([&, i](){
            for(auto&& p:traj_dict){
                auto& t = p.second;
                auto& id = p.first;
                if(id%hnsw_build_n_threads!=i){
                    continue;
                }
                // fmt::print("traj-{}, lens={}, #processed_traj={}\n", id, t.get_num_pnts(), cnt);
                sliding_window(t, sliding_window_length, sliding_window_step, sliding_window_interval, [&](const TrajView& tv){
                    int offset = std::distance(&t[0], (const Pnt3*)tv.p);
                    hnsw.addPoint(tv.p, SubtIdx(id, offset));
                    ++tot_cnt;
                });
            }
        });
    }
    for(int i=0;i<hnsw_build_n_threads;i++){
        thread_pool[i].join();
    }
    fmt::print("cnt={}\n", tot_cnt);
}

ObstacleDetector::~ObstacleDetector()
{

}


std::tuple<double, double, double, double> ObstacleDetector::get_densities(
    const TrajView& tv, 
    const TrajView& tv_succ
)
{
    std::vector<double> x = tv.to_flat_vec_spatial();
    std::vector<double> x_succ = tv_succ.to_flat_vec_spatial();

    double dq = get_density_q(x);
    double dq_succ = get_density_succ_q(x, x_succ);
    double dr = get_density_ref(x);
    double dr_succ = get_density_succ_ref(x, x_succ);
    return std::make_tuple(dq, dr, dq_succ, dr_succ);
}

//idx will be in ref trajectories
void ObstacleDetector::get_obstacle_one_subt(const SubtIdx& idx0)
{
    std::stack<SubtIdx> subtidx_stack;
    subtidx_stack.push(idx0);

    while(!subtidx_stack.empty()) {
        SubtIdx idx = subtidx_stack.top();
        subtidx_stack.pop();

        //checked before
        auto tmp_it = ref_is_candidates.find(idx);
        if(tmp_it != ref_is_candidates.end()){
            return ;
        }
        // TrajView tv_nr(ref_dict[idx.tid], idx.offset, sliding_window_len);
        // TrajView tv_nr_succ = try_get_trajview(ref_dict[idx.tid], idx.offset+1, sliding_window_len).value_or(tv_nr);
        auto [tv_nr, tv_nr_succ] = get_view_from_subtidx(*ref_dict_ptr, idx);

        auto [dq, dr, dq_comover, dr_comover] = get_densities(tv_nr, tv_nr_succ);
        // fmt::print("density={}, {}, {}, {}, q_checked.size()={}\n", dq, dr, dq_comover, dr_comover, q_checked.size());

        bool is_t_obstacle = obstacle_conditaion(dq, dr, dq_comover, dr_comover);
        ref_is_candidates[idx] = is_t_obstacle;

        if(is_t_obstacle){
            //set t as obstacle that belongs to q
            // Pnt3 lastp = tv_nr[tv_nr.get_num_pnts()-1];
            double ratio = (dr_comover/dr) / (dq_comover/dq);
            // fmt::print("ratio={}, densities={}, {}, {}, {}\n", ratio, dq, dr, dq_comover, dr_comover);
            nlohmann::json obsj = get_single_obstacle_json(dq, dr, dq_comover, dr_comover, tv_nr_succ);
            auto pidx = obst_subtidx_belongs_to.get(idx);
            if(obst_json_map.find(pidx)==obst_json_map.end()){ 
                obst_json_map[pidx] = json::array();
            }
            obst_json_map[pidx].push_back(obsj);

            int len = dim/2;

            hnsw_ref.forEachNeighbor0(idx, [&](const SubtIdx& qqidx, void* datap){
                double distqq = PntEuDist(datap, tv_nr.p, &len);
                //filter out far neighbors
                if(distqq > dist2_threshold){
                    return ;
                }
                merge_obstacle(idx, qqidx);
                if(q_checked.find(qqidx)!=q_checked.end()){
                    //already checked
                    return ;
                }
                subtidx_stack.push(qqidx);
            });
        }
    }
}

nlohmann::json ObstacleDetector::output_obstacle_json(const std::unordered_map<int, Traj>& query_dict)
{
    nlohmann::json j_ret, j_obst;
    
    MyTimer::pusht();
    // q_checked.clear();
    build_index_q(query_dict);
    obst_json_map.clear();

    sliding_window_trajdict(query_dict, [&](const SubtIdx& qidx, const TrajView& tv){
        if(q_checked.find(qidx) == q_checked.end()){
            q_checked.insert(qidx);
            // nlohmann::json obstacle_json;

            //get all neighbors of q
            // return a priority_queue<pair<double, SubtIdx> > 
            auto nr = hnsw_ref.searchKnn(tv.p, num_neighbors);
            for(;!nr.empty(); nr.pop()){
                auto nr_pair = nr.top();
                SubtIdx idx = nr_pair.second;
                if(nr_pair.first > dist2_threshold){
                    continue;
                }

                get_obstacle_one_subt(idx);
            }
        }
    });
    for(auto&& obst_json_pair:obst_json_map){ 
        if(obst_json_pair.second.size()>0){
            j_obst.emplace_back(std::move(obst_json_pair.second));
        }
    }
    double time_used = MyTimer::popt();
    j_ret["running_time"] = time_used;
    j_ret["obst"] = j_obst;
    return j_ret;
}



nlohmann::json ObstacleDetector::output_obstacle_json_with_index(const std::unordered_map<int, Traj>& query_dict, const std::string& filename)
{
    nlohmann::json j_ret, j_obst;
    
    MyTimer::pusht();
    // q_checked.clear();
    load_index_q(filename, query_dict);
    obst_json_map.clear();

    sliding_window_trajdict(query_dict, [&](const SubtIdx& qidx, const TrajView& tv){
        if(q_checked.find(qidx) == q_checked.end()){
            q_checked.insert(qidx);
            // nlohmann::json obstacle_json;

            //get all neighbors of q
            // return a priority_queue<pair<double, SubtIdx> > 
            auto nr = hnsw_ref.searchKnn(tv.p, num_neighbors);
            for(;!nr.empty(); nr.pop()){
                auto nr_pair = nr.top();
                SubtIdx idx = nr_pair.second;
                if(nr_pair.first > dist2_threshold){
                    continue;
                }

                get_obstacle_one_subt(idx);
            }
        }
    });
    for(auto&& obst_json_pair:obst_json_map){ 
        if(obst_json_pair.second.size()>0){
            j_obst.emplace_back(std::move(obst_json_pair.second));
        }
    }
    double time_used = MyTimer::popt();
    j_ret["running_time"] = time_used;
    j_ret["obst"] = j_obst;
    return j_ret;
}

void detect_obstacle(const std::unordered_map<int, Traj>& query_traj_dict)
{

}


void ObstacleDetector::save_index(const std::string& filename)
{
    std::string hnsw_filename = fmt::format("{}_hnsw.idx", filename);
    std::string kde_filename = fmt::format("{}_kde_sigma{}.idx", filename, sigma);
    std::string kde_succ_filename = fmt::format("{}_kde_succ_sigma{}.idx", filename, sigma);

    lsh_kde_ref.saveIndex(kde_filename);
    lsh_kde_succ_ref.saveIndex(kde_succ_filename);

    hnsw_ref.saveIndex(hnsw_filename);
}
void ObstacleDetector::save_index_q(const std::string& filename)
{
    std::string kde_filename = fmt::format("{}_kde_sigma{}.idx", filename, sigma);
    std::string kde_succ_filename = fmt::format("{}_kde_succ_sigma{}.idx", filename, sigma);

    lsh_kde_q.saveIndex(kde_filename);
    lsh_kde_succ_q.saveIndex(kde_succ_filename);
}
void ObstacleDetector::load_index(const std::string& filename, const std::unordered_map<int, Traj>& traj_dict)
{
    ref_dict_ptr = &traj_dict;

    std::string hnsw_filename = fmt::format("{}_hnsw.idx", filename);
    std::string kde_filename = fmt::format("{}_kde_sigma{}.idx", filename, sigma);
    std::string kde_succ_filename = fmt::format("{}_kde_succ_sigma{}.idx", filename, sigma);


    lsh_kde_ref.loadIndex(kde_filename);
    lsh_kde_succ_ref.loadIndex(kde_succ_filename);
    
    hnsw_ref.loadIndex(hnsw_filename, &hnsw_space, hnsw_max_objs);
}
void ObstacleDetector::load_index_q(const std::string& filename, const std::unordered_map<int, Traj>& traj_dict)
{
    queries_dict_ptr = &traj_dict;

    std::string kde_filename = fmt::format("{}_kde_sigma{}.idx", filename, sigma);
    std::string kde_succ_filename = fmt::format("{}_kde_succ_sigma{}.idx", filename, sigma);

    lsh_kde_q.loadIndex(kde_filename);
    lsh_kde_succ_q.loadIndex(kde_succ_filename);
}