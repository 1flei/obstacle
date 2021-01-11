#include "obstacle.h"


double ObsDetector::density(const SubtPriQueue& res_que)
{
    double ret = 0.;
    // fmt::print("            ---------------\n");
    for(const auto& p:res_que){
        // fmt::print("            ~{}\n", p);
        double dist2 = p.first;
        ret += exp(-dist2/sigma/sigma);
    }
    return ret;
}

double ObsDetector::density_prev(std::unordered_map<int, Traj>& traj_dict, const TrajView& tv_prev, 
    const SubtPriQueue& res_que, hnswlib::SpaceInterface<double>& distf)
{
    double ret = 0;
    for(const auto& p:res_que){
        int tid = p.second.tid;
        int offset = p.second.offset;

        TrajView tv_nr();
        Traj& t= traj_dict[tid];
        TrajView tv_p_prev = try_get_trajview(t, offset-1, swlen, swinterval).value_or(TrajView(t, offset, swlen));

        double dist2 = distf.get_dist_func()(tv_prev.p, tv_p_prev.p, distf.get_dist_func_param());
        ret += exp(-dist2/sigma/sigma);
    }
    return ret;
}

double ObsDetector::density_succ(std::unordered_map<int, Traj>& traj_dict, const TrajView& tv_succ, 
    const SubtPriQueue& res_que, hnswlib::SpaceInterface<double>& distf)
{
    // fmt::print("            ---------------\n");
    double ret = 0;
    for(const auto& p:res_que){
        int tid = p.second.tid;
        int offset = p.second.offset;

        TrajView tv_nr();
        Traj& t= traj_dict[tid];
        TrajView tv_p_succ = try_get_trajview(t, offset+1, swlen, swinterval).value_or(TrajView(t, offset, swlen));

        double dist2 = distf.get_dist_func()(tv_succ.p, tv_p_succ.p, distf.get_dist_func_param());
        ret += exp(-dist2/sigma/sigma);
        // fmt::print("            ~{}, {}\n", p, dist2);
    }
    return ret;
}

double ObsDetector::density_comover(std::unordered_map<int, Traj>& traj_dict, const TrajView& tv_succ, 
    const SubtPriQueue& res_que, hnswlib::SpaceInterface<double>& distf)
{
    // fmt::print("            ---------------\n");
    double ret = 0;
    for(const auto& p:res_que){
        int tid = p.second.tid;
        int offset = p.second.offset;

        TrajView tv_nr();
        Traj& t= traj_dict[tid];
        TrajView tv_p_succ = try_get_trajview(t, offset+1, swlen, swinterval).value_or(TrajView(t, offset, swlen));

        double dist2 = distf.get_dist_func()(tv_succ.p, tv_p_succ.p, distf.get_dist_func_param());
        double min_dist2 = std::max(dist2, p.first);
        ret += exp(-min_dist2/sigma/sigma);
        // fmt::print("            ~{}, {}\n", p, dist2);
    }
    return ret;
}

// double density(const SubtPriQueue& res_que, double sigma)
// {
//     double ret = 0.;
//     for(const auto& p:res_que){
//         double dist = p.first;
//         ret += exp(dist*dist/sigma/sigma);
//     }
//     return ret;
// }

double ObsDetector::get_score(double dq, double dr, double dq_succ, double dr_succ)
{
    return dr_succ / dr - dq_succ / dq;
}

double ObsDetector::get_score_fisher(double dq, double dr, double dq_comover, double dr_comover)
{
    double dq_ncom = dq-dq_comover;
    double dr_ncom = dr-dr_comover;
    return 0;
}

double ObsDetector::get_score_z(double dq, double dr, double dq_comover, double dr_comover)
{
    // double dq_ncom = dq-dq_comover;
    // double dr_ncom = dr-dr_comover;

    double p1_hat = dr_comover/dr;
    double p2_hat = dq_comover/dq;

    double p_hat = (dr*p1_hat + dq*p2_hat) / (dr+dq);
    double sigmad = sqrt(p_hat*(1-p_hat)*(1/dr + 1/dq));

    double z = (p1_hat-p2_hat) / sigmad;

    return z;
}


std::tuple<TrajView, TrajView> ObsDetector::get_view_from_subtidx(std::unordered_map<int, Traj>& traj_dict, const SubtIdx& idx)
{
    TrajView tv_nr(traj_dict[idx.tid], idx.offset, swlen);
    TrajView tv_nr_succ = try_get_trajview_clamp01(traj_dict[idx.tid], idx.offset+swstep, swlen, swinterval).value_or(tv_nr);
    return std::make_tuple(tv_nr, tv_nr_succ);
}


void ObsDetector::filter_res_que(
    std::unordered_map<int, Traj>& traj_dict, 
    SubtPriQueue& res_que
)
{
    for(auto it=res_que.begin(); it!=res_que.end();){
        const auto& p = *it;
        int tid = p.second.tid;
        int offset = p.second.offset;

        Traj& t= traj_dict[tid];
        auto maybe_tv_succ = try_get_trajview(t, offset+1, swlen, swinterval);
        if(!maybe_tv_succ){
            //consider only those trajectories with next subt
            it = res_que.erase(it);
        } else{
            ++it;
        }
    }
}

std::tuple<double, double> ObsDetector::get_densities(
    std::unordered_map<int, Traj>& traj_dict, 
    const TrajView& tv_succ, 
    const SubtPriQueue& res_que, 
    hnswlib::SpaceInterface<double>& distf, 
    double dist_sigma
)
{
    double d_cur = 0.;
    double d_comover = 0.;
    for(const auto& p:res_que){
        int tid = p.second.tid;
        int offset = p.second.offset;

        Traj& t= traj_dict[tid];
        auto maybe_tv_succ = try_get_trajview(t, offset+1, swlen, swinterval);
        if(!maybe_tv_succ){
            //consider only those trajectories with next subt
            continue;
        }
        const TrajView& tv_p_succ = *maybe_tv_succ;

        double dist2 = p.first;
        double dist2_succ = distf.get_dist_func()(tv_succ.p, tv_p_succ.p, distf.get_dist_func_param());
        double min_dist2 = std::max(dist2, dist2_succ);
        d_cur += exp(-dist2/dist_sigma/dist_sigma);
        d_comover += exp(-min_dist2/dist_sigma/dist_sigma);
        // fmt::print("            ~{}, {}\n", p, dist2);
    }
    return std::make_tuple(d_cur, d_comover);
}

std::tuple<double, double, double, double> ObsDetector::get_densities(
    const TrajView& tv_nr, 
    const TrajView& tv_nr_succ, 
    SubtPriQueue& nqr, 
    SubtPriQueue& nrr, 
    std::unordered_map<int, Traj>& ref_dict, 
    std::unordered_map<int, Traj>& queries_dict, 
    hnswlib::SpaceInterface<double>& dtwspace
)
{
    // auto it = nrr.begin();
    // std::advance(it, nrr.size()/2);
    // double median_dist_in_nr = sqrt(it->first);

    // double min_dist = sqrt(nrr.minimum().first);
    // double max_dist = sqrt(nrr.maximum().first);
    // fmt::print("median_dist_in_nr={}, min_dist={}, max_dist={}\n", median_dist_in_nr, min_dist, max_dist);
    // double dist_sigma = std::min(max_dist, SIGMA);
    double dist_sigma = sigma;


    auto [dq, dq_comover] = get_densities(queries_dict, tv_nr_succ, nqr, dtwspace, dist_sigma);
    auto [dr, dr_comover] = get_densities(ref_dict, tv_nr_succ, nrr, dtwspace, dist_sigma);
    return std::make_tuple(dq, dr, dq_comover, dr_comover);
}

nlohmann::json ObsDetector::get_single_obstacle_json(
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
    const SubtPriQueue& nrr
)
{
    nlohmann::json j;
    // fmt::print("tv_nr_succ={}, its.num_pnts()={}\n", tv_nr_succ, tv_nr_succ.get_num_pnts());
    Pnt3 lastp = tv_nr_succ[tv_nr_succ.get_num_pnts()-1];
    std::vector<double> tmp{lastp.t, lastp.x, lastp.y};

    j["obs_pnts"] = tmp;

    // fmt::print("1!!!!----{}\n", j.dump(4));

    if(verbosity>=1) {
        j["res_p_dist"] = res_p_dist;
        j["res_p_tid"] = nr_idx.tid;
        j["res_p_offset"] = nr_idx.offset;

        std::vector<double> densities{dq, dr, dq_comover, dr_comover};
        j["densities"] = densities;
    }
    // fmt::print("2!!!!----{}\n", j.dump(4));

    if(verbosity>=2){
        auto [tv_q, tv_q_succ] = get_view_from_subtidx(queries_dict, qidx);
        j["q"] = tv_q.to_flat_vec();
        j["q_succ"] = tv_q_succ.to_flat_vec();
        j["t"] = tv_nr.to_flat_vec();
        j["t_succ"] = tv_nr_succ.to_flat_vec();
    }
    // fmt::print("3!!!!----{}\n", j.dump(4));

    if(verbosity>=3) {
        for(auto& p:nqr){
            if(p.first > dist2_threshold){
                continue;
            }
            auto [nqt_tv, nqt_tv_succ] = get_view_from_subtidx(queries_dict, p.second);
            j["Nqt"].push_back(nqt_tv.to_flat_vec());
            j["Nqt_succ"].push_back(nqt_tv_succ.to_flat_vec());
        }
        for(auto& p:nrr){
            if(p.first > dist2_threshold){
                continue;
            }
            auto [nrt_tv, nrt_tv_succ] = get_view_from_subtidx(ref_dict, p.second);
            j["Nrt"].push_back(nrt_tv.to_flat_vec());
            j["Nrt_succ"].push_back(nrt_tv_succ.to_flat_vec());
        }
    }
    // fmt::print("4!!!!----{}\n", j.dump(4));
    return j;
}

void ObsDetector::get_obstacle_one_subt(
    const SubtIdx& qidx, 
    const TrajView& tv,  
    nlohmann::json& obstacle_json
)
{
    q_checked.insert(qidx);

    //get all neighbors of q
    auto nr = searchKNN(hnsw, tv.p, numNeighbors);

    for(const auto& nr_pair:nr){
        double distq = nr_pair.first;
        SubtIdx idx = nr_pair.second;
        
        auto tmp_it = ref_is_candidates.find(idx);
        if(tmp_it != ref_is_candidates.end()){
            continue;
        }
        // TrajView tv_nr(ref_dict[idx.tid], idx.offset, sliding_window_len);
        // TrajView tv_nr_succ = try_get_trajview(ref_dict[idx.tid], idx.offset+1, sliding_window_len).value_or(tv_nr);
        auto [tv_nr, tv_nr_succ] = get_view_from_subtidx(ref_dict, idx);
        if(nr_pair.first > dist2_threshold){
            // fmt::print("dist={}\n", nr_pair.first);
            continue;
        }

        auto nqr = searchKNN(hnsw_q, tv_nr.p, numNeighbors);
        nqr.emplace_k(distq, qidx, numNeighbors);

        auto nrr = [&](){
            auto ref_neighbor_cache_it = ref_neighbor_cache.find(idx);
            if(ref_neighbor_cache_it!=ref_neighbor_cache.end()){
                SubtPriQueue que;
                for(auto &p:ref_neighbor_cache_it->second) {
                    que.emplace_k(p.first, p.second, numNeighbors);
                }
                return que;
            } else{
                auto que = searchKNN(hnsw, tv_nr.p, numNeighbors);
                for(auto& p:que){
                    ref_neighbor_cache[idx].push_back(p);
                }
                return que;
            }
        }();
        filter_res_que(queries_dict, nqr);
        filter_res_que(ref_dict, nrr);

        // auto nqr_succ = searchKNN(hnsw_q, tv_nr_succ.p, numNeighbors);
        // auto nrr_succ = searchKNN(hnsw, tv_nr_succ.p, numNeighbors);

        auto [dq, dr, dq_comover, dr_comover] = get_densities(tv_nr, tv_nr_succ, nqr, nrr, ref_dict, queries_dict, dtwspace);
        // fmt::print("density={}, {}, {}, {}, q_checked.size()={}\n", dq, dr, dq_comover, dr_comover, q_checked.size());

        bool is_t_obstacle = obstacle_conditaion(dq, dr, dq_comover, dr_comover);
        ref_is_candidates[idx] = is_t_obstacle;

        //opt-3
        for(auto& nrr_pair: nrr){
            double distr = nrr_pair.first;
            if(distr < 0.05*sigma){
                SubtIdx ridx = nrr_pair.second;
                //then r will have the same type of r
                ref_is_candidates[ridx] = is_t_obstacle;
                if(is_t_obstacle){
                    auto [r_tv_nr, r_tv_nr_succ] = get_view_from_subtidx(ref_dict, ridx);
                    nlohmann::json obsj = get_single_obstacle_json(r_tv_nr, r_tv_nr_succ, distq, ridx, qidx, dq, dr, dq_comover, dr_comover, nqr, nrr);
                    obstacle_json.push_back(obsj);
                }
            }
        }

        if(is_t_obstacle){
            //set t as obstacle that belongs to q
            // Pnt3 lastp = tv_nr[tv_nr.get_num_pnts()-1];
            nlohmann::json obsj = get_single_obstacle_json(tv_nr, tv_nr_succ, distq, idx, qidx, dq, dr, dq_comover, dr_comover, nqr, nrr);
            obstacle_json.push_back(obsj);
            auto nq = searchKNN(hnsw_q, tv.p, numNeighbors);
            // auto& nq = nqr;

            for(auto& nqr_pair:nq) {
                double distqq = nqr_pair.first;
                SubtIdx qqidx = nqr_pair.second;
                if(q_checked.find(qqidx)!=q_checked.end()){
                    //already checked
                    continue ;
                }
                if(distqq < 0.05*sigma){
                    continue ;
                }
                auto [qq_tv, qq_tv_succ] = get_view_from_subtidx(queries_dict, qqidx);
                get_obstacle_one_subt(qqidx, qq_tv, obstacle_json);
            }
        }
    }
}

void ObsDetector::get_obstacle_one_subt_naive(
    const SubtIdx& qidx, 
    const TrajView& tv,  
    nlohmann::json& obstacle_json
)
{
    q_checked.insert(qidx);

    //get all neighbors of q
    auto nr = searchKNN(hnsw, tv.p, numNeighbors);

    for(const auto& nr_pair:nr){
        double distq = nr_pair.first;
        SubtIdx idx = nr_pair.second;
        
        // TrajView tv_nr(ref_dict[idx.tid], idx.offset, sliding_window_len);
        // TrajView tv_nr_succ = try_get_trajview(ref_dict[idx.tid], idx.offset+1, sliding_window_len).value_or(tv_nr);
        auto [tv_nr, tv_nr_succ] = get_view_from_subtidx(ref_dict, idx);
        if(nr_pair.first > dist2_threshold){
            // fmt::print("dist={}\n", nr_pair.first);
            continue;
        }

        auto nqr = searchKNN(hnsw_q, tv_nr.p, numNeighbors);
        nqr.emplace_k(distq, qidx, numNeighbors);
        auto nrr = searchKNN(hnsw, tv_nr.p, numNeighbors);
        filter_res_que(queries_dict, nqr);
        filter_res_que(ref_dict, nrr);

        // auto nqr_succ = searchKNN(hnsw_q, tv_nr_succ.p, numNeighbors);
        // auto nrr_succ = searchKNN(hnsw, tv_nr_succ.p, numNeighbors);

        auto [dq, dr, dq_comover, dr_comover] = get_densities(tv_nr, tv_nr_succ, nqr, nrr, ref_dict, queries_dict, dtwspace);
        // fmt::print("density={}, {}, {}, {}, q_checked.size()={}\n", dq, dr, dq_comover, dr_comover, q_checked.size());

        bool is_t_obstacle = obstacle_conditaion(dq, dr, dq_comover, dr_comover);
        if(is_t_obstacle){
            //set t as obstacle that belongs to q
            // Pnt3 lastp = tv_nr[tv_nr.get_num_pnts()-1];
            nlohmann::json obsj = get_single_obstacle_json(tv_nr, tv_nr_succ, distq, idx, qidx, dq, dr, dq_comover, dr_comover, nqr, nrr);
            obstacle_json.push_back(obsj);
    
            auto nq = searchKNN(hnsw_q, tv.p, numNeighbors);
            for(auto& nqr_pair:nq) {
                double distqq = nqr_pair.first;
                SubtIdx qqidx = nqr_pair.second;
                if(q_checked.find(qqidx)!=q_checked.end()){
                    //already checked
                    continue ;
                }
                auto [qq_tv, qq_tv_succ] = get_view_from_subtidx(queries_dict, qqidx);
                get_obstacle_one_subt(qqidx, qq_tv, obstacle_json);
            }
        }
    }
}

nlohmann::json ObsDetector::output_obstacle_json()
{
    nlohmann::json j_ret, j_obst;
    
    MyTimer::pusht();
    sliding_window_trajdict(queries_dict, [&](const SubtIdx& qidx, const TrajView& tv){
        if(q_checked.find(qidx) == q_checked.end()){
            nlohmann::json obstacle_json;
            get_obstacle_one_subt(qidx, tv, obstacle_json);
            
            if(obstacle_json.size() > 0){
                // fmt::print("qtv={}, #obs_pnts={}\n", tv.to_flat_vec(), obstacle_json.size());
                j_obst.emplace_back(std::move(obstacle_json));
            }
        }
    });
    double time_used = MyTimer::popt();
    j_ret["running_time"] = time_used;
    j_ret["obst"] = j_obst;
    return j_ret;
}

nlohmann::json ObsDetector::output_obstacle_json_naive()
{
    nlohmann::json j_ret, j_obst;
    
    MyTimer::pusht();
    sliding_window_trajdict(queries_dict, [&](const SubtIdx& qidx, const TrajView& tv){
        if(q_checked.find(qidx) == q_checked.end()){
            nlohmann::json obstacle_json;
            get_obstacle_one_subt_naive(qidx, tv, obstacle_json);
            
            if(obstacle_json.size() > 0){
                // fmt::print("qtv={}, #obs_pnts={}\n", tv.to_flat_vec(), obstacle_json.size());
                j_obst.emplace_back(std::move(obstacle_json));
            }
        }
    });
    double time_used = MyTimer::popt();
    j_ret["running_time"] = time_used;
    j_ret["obst"] = j_obst;
    return j_ret;
}

nlohmann::json ObsDetector::output_subtraj_density_json_bbox(const BBox& bb)
{
    int querycnt = 0;

    nlohmann::json j_ret, j_obst;
    sliding_window_trajdict(queries_dict, [&](const SubtIdx& qidx, const TrajView& tv){
        // fmt::print("query-{}, qidx={}, trajs={}\n", querycnt, qidx, tv);
        auto& p = tv[tv.get_num_pnts()-1];
        if(p.x > bb.maxlng || p.x < bb.minlng || p.y > bb.maxlat || p.y < bb.minlat || p.t < bb.mint || p.t> bb.maxt) {
            return ;
        }

        fmt::print("p={}, tv={}, tv.num_pnts()={}\n", p, tv, tv.get_num_pnts());

        auto nr = searchKNN(hnsw, tv.p, numNeighbors);

        for(const auto& nr_pair:nr){
            double distq = nr_pair.first;
            SubtIdx idx = nr_pair.second;
            if(ref_is_candidates.find(idx)!=ref_is_candidates.end()){
                continue;
            }
            // TrajView tv_nr(ref_dict[idx.tid], idx.offset, sliding_window_len);
            // TrajView tv_nr_succ = try_get_trajview(ref_dict[idx.tid], idx.offset+1, sliding_window_len).value_or(tv_nr);
            auto [tv_nr, tv_nr_succ] = get_view_from_subtidx(ref_dict, idx);
            if(nr_pair.first > dist2_threshold){
                continue;
            }


            auto nqr = searchKNN(hnsw_q, tv_nr.p, numNeighbors);
            auto nrr = searchKNN(hnsw, tv_nr.p, numNeighbors);
            nqr.emplace_k(distq, qidx, numNeighbors);
            // auto nqr_succ = searchKNN(hnsw_q, tv_nr_succ.p, numNeighbors);
            // auto nrr_succ = searchKNN(hnsw, tv_nr_succ.p, numNeighbors);
            filter_res_que(queries_dict, nqr);
            filter_res_que(ref_dict, nrr);


            auto [dq, dr, dq_comover, dr_comover] = get_densities(tv_nr, tv_nr_succ, nqr, nrr, ref_dict, queries_dict, dtwspace);

            double score_z = get_score_z(dq, dr, dq_comover, dr_comover);
            fmt::print("density={}, {}, {}, {}; z-score={}\n", dq, dr, dq_comover, dr_comover, score_z);

            nlohmann::json j = get_single_obstacle_json(tv_nr, tv_nr_succ, distq, idx, qidx, dq, dr, dq_comover, dr_comover, nqr, nrr);

            j_obst.push_back(j);

            ref_is_candidates[idx] = obstacle_conditaion(dq, dr, dq_comover, dr_comover);
        }
        querycnt++;
    });
    j_ret["obst"] = j_obst;
    return j_ret;
}