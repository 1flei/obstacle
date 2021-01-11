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


//detect anomaly trajectory from traj_dict
nlohmann::json kde_anomaly(
    std::unordered_map<int, Traj>& traj_dict, 
    ObstacleDetector& obstacle_detector, 
    double threshold)
{
    nlohmann::json j_ret;
    //for each data object, test the true density and the estimated density via lsh_sample
    for(auto&& p:traj_dict){
        auto& t = p.second;
        auto& id = p.first;
        obstacle_detector.for_consecutive_subtraj_spatial_vec(t, [&](const std::vector<double>& q, const std::vector<double>& q_succ){
            double est_kde = obstacle_detector.get_density_ref(q);
            if(est_kde < threshold) {
                //an anomaly
                nlohmann::json j;
                j["flat_spatial_array"] = q;
                j["est_kde"] = est_kde;
                // fmt::print("id={}, vec={}, est_kde={}\n", id, q, est_kde);

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

    double density_threshold = 10;

    std::vector<double> bbox_args;
    
    int maxObjs = 6000000;
    size_t M = 16;
    size_t ef_construction = 200;
    int numNeighbors = 8;
    int n_thread_building_index = 16;

    int lsh_density_nrepeat = 64;


    srand(666);
    options_description desc("Allowed options");
    desc.add_options()
        ("help,h", "help message")

        //parameters of HNSW, might use default value
        ("hnsw_max_objs", value(&maxObjs), "the maximum number of objects in HNSW")
        ("M", value(&M), "M for HNSW")
        ("ef_construction", value(&ef_construction), "ef_construction for HNSW")


        //model paramters
        ("sigma", value(&sigma), "sigma for KDE")
        ("significance", value(&significance), "significance level of score")
        ("min_support", value(&min_support), "minimum support to be an obstacle")
        ("n_thread_building_index", value(&n_thread_building_index), "number of threads used for building index")
        ("threshold", value(&density_threshold), "density threshold")

        ("numNeighbors,K", value(&numNeighbors), "number of neighbors considered")
        ("verbolity_level,V", value(&verbosity_level), "verbosity_level")
        ("dataColName,D", value(&dataCollectionName), "data collection name")
        ("queryColName,Q", value(&queryCollectionName), "query collection name")
        ("outputFileName,O", value(&output_filename), "output file name (usually xxx.json)")

        ("time_interval", value(&time_interval), "time interval used")
        ("sliding_window_len", value(&sliding_window_len), "the length of sliding windows")
        ("sliding_window_step", value(&sliding_window_len), "the step of sliding windows")
        
        ("lsh_density_nrepeat", value(&lsh_density_nrepeat), "time interval used")

        ("dataIndexFileName", value(&dataIndexFIle), "indexing file name")
        ("queryIndexFileName", value(&queryIndexFile), "query indexing file name")

        ("is_visualizing_subt", "whether to use the function of visualizing sub-traj")
        ("naive_version", "version version that is used for the experiment only")
        ("bbox", value(&bbox_args)->multitoken(), "bounding box used")

        ("force_build_index", "force rebuild index")
    ;
    variables_map vm;
    store(parse_command_line(argc, argv, desc), vm);
    notify(vm);  

    if (vm.count("help")) {
        std::cout << desc << std::endl;
        return 1;
    }

    std::string index_filename = fmt::format("{}_kde_sigma{}.idx", dataCollectionName, sigma);

    auto ref_dict = read_flat_trajs(conn, dataCollectionName);
    auto queries_dict = read_flat_trajs(conn, queryCollectionName);
    fmt::print("ref_dict.size()={}, queries_dict.size()={}\n", ref_dict.size(), queries_dict.size());   

    ObstacleDetector obstacle_detector(lsh_density_nrepeat, sliding_window_len, sliding_window_step, time_interval, 
        sigma, min_support, significance, maxObjs, M, ef_construction, n_thread_building_index);
    fmt::print("building index\n");
    auto indexing_time = MyTimer::funcTime([&](){
        if(vm.count("force_build_index")){
            fmt::print("force to build from databases\n");

            obstacle_detector.build_lsh_density_ref(ref_dict);
        } else{
            try {
                //for test purpose
                obstacle_detector.get_lsh_kde_ref().loadIndex(index_filename);
            } catch(const std::exception& e) {
                fmt::print("exception when reading index: {}\n", e.what());
                fmt::print("so build from databases\n");

                obstacle_detector.build_lsh_density_ref(ref_dict);
                obstacle_detector.get_lsh_kde_ref().saveIndex(index_filename);
            }
        }
    });

    fmt::print("testing\n");
    std::ofstream outf(output_filename);
    auto j_ret = kde_anomaly(queries_dict, obstacle_detector, density_threshold);
    outf << j_ret.dump(4);

    return 0;
}