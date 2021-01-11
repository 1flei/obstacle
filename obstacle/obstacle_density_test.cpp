#include "boost/program_options.hpp"
#include "../mongodb/driver.hpp"
#include <iostream>
#include "../my_timer.h"
#include "../lsh_density/lsh_density.h"
#include "../lsh_density/lsh_density_linear.h"
#include "json.hpp"
#include "obstacle_detector.h"
#include <tuple>

using namespace boost::program_options;


std::tuple<double, double> brute_force_kde_succ(std::unordered_map<int, Traj>& traj_dict, 
    const std::vector<double>& vec_q, 
    const std::vector<double>& vec_q_succ, 
    ObstacleDetector& obstacle_detector)
{
    double kde = 0;
    double kde_succ = 0;
    int nsubt = 0;
    double sigma = obstacle_detector.get_sigma();
    int len = obstacle_detector.get_sliding_window_length();
    for(auto&& pdata: traj_dict) {
        auto& t = pdata.second;
        auto& id = pdata.first;

        obstacle_detector.for_consecutive_subtraj_spatial_vec(t, [&](const std::vector<double>& x, const std::vector<double>& x_succ){
            double dist2x = calc_l2_sqr<double>(len*2, &vec_q[0], &x[0]);
            double dist2x_succ = calc_l2_sqr<double>(len*2, &vec_q_succ[0], &x_succ[0]);
            kde += exp(-dist2x * 0.5 / sigma/sigma);
            kde_succ += exp(-dist2x * 0.5 / sigma/sigma)*exp(-dist2x_succ * 0.5 / sigma/sigma);
            nsubt++;
        });
    }
    return std::make_tuple(kde, kde_succ);
}

nlohmann::json test_density(std::unordered_map<int, Traj>& traj_dict_q, std::unordered_map<int, Traj>& traj_dict, ObstacleDetector& obstacle_detector)
{
    nlohmann::json j_ret;
    //for each data object, test the true density and the estimated density via lsh_sample
    for(auto&& p:traj_dict_q){
        auto& t = p.second;
        auto& id = p.first;
        obstacle_detector.for_consecutive_subtraj_spatial_vec(t, [&](const std::vector<double>& q, const std::vector<double>& q_succ){
            //test only 1%
            if(rand()%10000<=100){
                nlohmann::json j;
                double est_kde = obstacle_detector.get_density_ref(q);
                double est_kde_succ = obstacle_detector.get_density_succ_ref(q, q_succ);
                j["flat_spatial_array"] = q;

                auto[kde, kde_succ] = brute_force_kde_succ(traj_dict, q, q_succ, obstacle_detector);
                j["est_kde"] = est_kde;
                j["brute_force_kde"] = kde;
                j["est_kde_succ"] = est_kde_succ;
                j["brute_force_kde_succ"] = kde_succ;

                fmt::print("id={}, vec={}, est_kde={}, kde={}, est_kde_succ={}, kde_succ={}\n", id, q, est_kde, kde, est_kde_succ, kde_succ);

                j_ret.push_back(j);
            }
        });
    }
    return j_ret;
}

int main(int argc, char **argv)
{
    mongocxx::instance inst{};
    mongocxx::client conn{mongocxx::uri{}};
    std::string dataCollectionName = "sgtaxi_ref";
    std::string queryCollectionName = "sgtaxi_afternoon_q";

    double time_interval = 30;  //30 seconds
    int sliding_window_len = 6;
    int sliding_window_step = 1;

    double sigma = 0.001;
    double significance = 1.645;
    double min_support = 2.;

    std::string output_filename = "obst_test.json";

    std::string dataIndexFIle;
    std::string queryIndexFile;

    int verbosity_level = 3;

    int k = 8;
    int l = 3;
    int nrepeat = 64;
    int count_min_nrepeat = 5;
    int bottom_k_sketch_size = 32;
    double r = 5.12;

    int ht_size = 4194301;

    std::vector<double> bbox_args;
    
    int maxObjs = 6000000;
    size_t M = 16;
    size_t ef_construction = 200;
    int numNeighbors = 8;
    int n_thread_building_index = 16;


    srand(666);
    options_description desc("Allowed options");
    desc.add_options()
        ("help,h", "help message")

        ("k,k", value(&k), "k for LSH")
        ("l,l", value(&l), "l for LSH")
        ("nrepeat,t", value(&nrepeat), "# repeated for the LSH_KDE")
        ("count_min_nrepeat,m", value(&count_min_nrepeat), "# repeated for the count-min sketch")
        ("ht_size,H", value(&ht_size), "size of each individual hash table")

        //model paramters
        ("sigma", value(&sigma), "sigma for KDE")
        ("significance", value(&significance), "significance level of score")
        ("min_support", value(&min_support), "minimum support to be an obstacle")

        ("verbolity_level,V", value(&verbosity_level), "verbosity_level")
        ("dataColName,D", value(&dataCollectionName), "data collection name")
        ("queryColName,Q", value(&queryCollectionName), "query collection name")
        ("outputFileName,O", value(&output_filename), "output file name (usually xxx.json)")

        ("time_interval", value(&time_interval), "time interval used")
        ("sliding_window_len", value(&sliding_window_len), "the length of sliding windows")
        ("sliding_window_step", value(&sliding_window_len), "the step of sliding windows")

        ("dataIndexFileName", value(&dataIndexFIle), "indexing file name")
        ("queryIndexFileName", value(&queryIndexFile), "query indexing file name")

        ("is_visualizing_subt", "whether to use the function of visualizing sub-traj")
        ("bbox", value(&bbox_args)->multitoken(), "bounding box used")

        ("force_build_index", "force rebuild index")

        //parameters of HNSW, might use default value
        ("maxObjs", value(&maxObjs), "the maximum number of objects in HNSW")
        ("M", value(&M), "M for HNSW")
        ("ef_construction", value(&ef_construction), "ef_construction for HNSW")
    ;
    variables_map vm;
    store(parse_command_line(argc, argv, desc), vm);
    notify(vm);  

    if (vm.count("help")) {
        std::cout << desc << std::endl;
        return 1;
    }

    if(dataIndexFIle==""){
        dataIndexFIle = fmt::format("{}.lshidx", dataCollectionName);
    }
    if(queryIndexFile==""){
        queryIndexFile = fmt::format("{}.lshidx", queryCollectionName);
    }


    auto ref_dict = read_flat_trajs(conn, dataCollectionName);
    auto queries_dict = read_flat_trajs(conn, queryCollectionName);
    fmt::print("ref_dict.size()={}, queries_dict.size()={}\n", ref_dict.size(), queries_dict.size());   

    ObstacleDetector obstacle_detector(nrepeat, sliding_window_len, sliding_window_step, time_interval, sigma, significance, min_support, maxObjs, M, ef_construction, n_thread_building_index);
    fmt::print("building index\n");
    obstacle_detector.build_lsh_density_ref(ref_dict);

    fmt::print("testing\n");
    std::ofstream outf(output_filename);
    auto j_ret = test_density(queries_dict, ref_dict, obstacle_detector);
    outf << j_ret.dump(4);

    return 0;
}