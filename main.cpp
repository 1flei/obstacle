#include "index/index.h"
#include "mongodb/driver.hpp"
#include <iostream>
// #include "obstacle/obstacle.h"
#include "boost/program_options.hpp"
#include "my_timer.h"
#include "obstacle/obstacle_detector.h"

using namespace boost::program_options;

int main(int argc, char **argv)
{
    //read data/query
    //data/data_index should be serialized before

    mongocxx::instance inst{};
    mongocxx::client conn{mongocxx::uri{}};

    int maxObjs = 6000000;
    size_t M = 16;
    size_t ef_construction = 200;
    int numNeighbors = 32;
    int n_thread_building_index = 4;

    int verbosity_level = 3;
    std::string dataCollectionName = "sgtaxi_ref";
    std::string queryCollectionName = "sgtaxi_afternoon_q";

    double time_interval = 30;  //30 seconds
    int sliding_window_len = 6;
    int sliding_window_step = 1;

    double sigma = 0.5;
    double significance = 1.645;
    double min_support = 2.;

    Pnt3DTWNormalizedSpace dtwspace(sliding_window_len);

    std::string output_filename = "obst_test.json";

    std::string dataIndexFIle;
    std::string queryIndexFile;

    int lsh_density_nrepeat = 64;


    std::vector<double> bbox_args;

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


    auto ref_dict = read_flat_trajs(conn, dataCollectionName);
    auto queries_dict = read_flat_trajs(conn, queryCollectionName);
    fmt::print("ref_dict.size()={}, queries_dict.size()={}\n", ref_dict.size(), queries_dict.size());   

    std::ofstream outf(output_filename);
    
    ObstacleDetector obstacle_detector(lsh_density_nrepeat, sliding_window_len, sliding_window_step, time_interval, 
        sigma, min_support, significance, maxObjs, M, ef_construction, n_thread_building_index);
    auto indexing_time = MyTimer::funcTime([&](){
        if(vm.count("force_build_index")){
            fmt::print("force to build from databases\n");

            obstacle_detector.build_index_ref(ref_dict);
        } else{
            try {
                //for test purpose
                obstacle_detector.load_index(dataCollectionName, ref_dict);
            } catch(const std::exception& e) {
                fmt::print("exception when reading index: {}\n", e.what());
                fmt::print("so build from databases\n");

                obstacle_detector.build_index_ref(ref_dict);
                obstacle_detector.save_index(dataCollectionName);
            }
        }
    });

    nlohmann::json j_ret;
    if(vm.count("force_build_index")){
        j_ret = obstacle_detector.output_obstacle_json(queries_dict);
    } else {
        try {
            j_ret = obstacle_detector.output_obstacle_json_with_index(queries_dict, queryCollectionName);
        } catch(const std::exception& e) {
            fmt::print("exception when reading index_q: {}\n", e.what());
            fmt::print("so build from databases\n");

            j_ret = obstacle_detector.output_obstacle_json(queries_dict);
            obstacle_detector.save_index_q(queryCollectionName);
        }
    }
    j_ret["indexing_time"] = indexing_time.count();
    j_ret["num_ref_pnts"] = obstacle_detector.get_num_ref_pnts();
    j_ret["num_q_pnts"] = obstacle_detector.get_num_q_pnts();
    
    outf << j_ret.dump(4);

    return 0;
}