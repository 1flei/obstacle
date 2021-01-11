#include "../lsh_density/lsh_density_bks.h"
#include "boost/program_options.hpp"
#include "../mongodb/driver.hpp"
#include <iostream>
#include "../my_timer.h"
#include "index.h"
#include "json.hpp"

using namespace boost::program_options;


double brute_force_kde(std::unordered_map<int, Traj>& traj_dict, const TrajView& tv_q, double sigma, int len, int step, double interval)
{
    std::vector<double> vec_q = tv_q.to_flat_vec_spatial();
    double kde = 0;
    for(auto&& pdata: traj_dict) {
        auto& t = pdata.second;
        auto& id = pdata.first;
        sliding_window(t, len, step, interval, [&](const TrajView& tv_data){
            std::vector<double> vec_di = tv_data.to_flat_vec_spatial();
            double dist2 = calc_l2_sqr<double>(len*2, &vec_q[0], &vec_di[0]);
            kde += exp(-dist2 * 0.5 / sigma/sigma);
        });
    }
    return kde;
}

nlohmann::json test_density(std::unordered_map<int, Traj>& traj_dict_q, std::unordered_map<int, Traj>& traj_dict, LSH_KDE& index, double sigma, 
        int len, int step, double interval)
{
    nlohmann::json j_ret;
    //for each data object, test the true density and the estimated density via lsh_sample
    for(auto&& p:traj_dict_q){
        auto& t = p.second;
        auto& id = p.first;
        sliding_window(t, len, step, interval, [&](const TrajView& tv){
            //test only 1%
            if(rand()%10000<=50){
                nlohmann::json j;
                std::vector<double> flat_spatial_vec = tv.to_flat_vec_spatial();
                double est_kde = index.get_estimated_density(flat_spatial_vec);
                j["flat_spatial_array"] = flat_spatial_vec;

                double kde = brute_force_kde(traj_dict, tv, sigma, len, step, interval);
                j["lsh_cnt"] = est_kde;
                j["brute_force_kde"] = kde;

                fmt::print("id={}, vec={}, lsh_cnt={}, kde={}\n", id, flat_spatial_vec, est_kde, kde);

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

    LSH_KDE index(sliding_window_len*2, k, l, sigma*r, nrepeat, ht_size, count_min_nrepeat);
    LSH_KDE index_q(sliding_window_len*2, k, l, sigma*r, nrepeat, ht_size, count_min_nrepeat);

    //then build lsh_sample for it

    auto indexing_time = MyTimer::funcTime([&](){
        if(vm.count("force_build_index")){
            fmt::print("force to build from databases\n");

            build_index(ref_dict, sliding_window_len, sliding_window_step, time_interval, index);
            index.saveIndex(dataIndexFIle);
        } else{
            try {
                index.loadIndex(dataIndexFIle);
            } catch(const std::exception& e) {
                fmt::print("exception when reading index: {}\n", e.what());
                fmt::print("so build from databases\n");

                build_index(ref_dict, sliding_window_len, sliding_window_step, time_interval, index);
                index.saveIndex(dataIndexFIle);
            }
        }
    });

    auto query_indexing_time = MyTimer::funcTime([&](){
        if(vm.count("force_build_index")){
            fmt::print("force to build from databases\n");

            build_index(queries_dict, sliding_window_len, sliding_window_step, time_interval, index_q);
            index_q.saveIndex(queryIndexFile);
        } else{
            try {
                //for test purpose
                index_q.loadIndex(queryIndexFile);
            } catch(const std::exception& e) {
                fmt::print("exception when reading index: {}\n", e.what());
                fmt::print("so build from databases\n");

                build_index(queries_dict, sliding_window_len, sliding_window_step, time_interval, index_q);
                index_q.saveIndex(queryIndexFile);
            }
        }
    });
    //build query index
    // build_index(queries_dict, hnsw_q);

    std::ofstream outf(output_filename);
    auto j_ret = test_density(queries_dict, ref_dict, index, sigma, sliding_window_len, sliding_window_step, time_interval);
    outf << j_ret.dump(4);

    return 0;
}