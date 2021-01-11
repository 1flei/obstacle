#include <string>
#include "../mongodb/driver.hpp"
#include <iostream>
#include "../my_timer.h"
#include "index.h"
#include "json.hpp"
#include "dataset_io.h"

void gen_datasets(mongocxx::client& conn, 
    const std::string& database_collection_name, 
    const std::string& output_filename,
    int len, int step, double interval, 
    double sample_rate=1.001)
{
    std::unordered_map<int, Traj> traj_dict = read_flat_trajs(conn, database_collection_name);
    std::vector<std::vector<double> > datasets; 
    for(auto&& p:traj_dict){
        auto& t = p.second;
        auto& id = p.first;
        sliding_window(t, len, step, interval, [&](const TrajView& tv){
            //test only 1%
            if((rand()%10000)<=10000*sample_rate){
                std::vector<double> flat_spatial_vec = tv.to_flat_vec_spatial();
                datasets.emplace_back(std::move(flat_spatial_vec));
            }
        });
    }

    write_dataset(output_filename, datasets);
}

//generate dataset for sgtaxi and the vessel dataset


int main()
{
    mongocxx::instance inst{};
    mongocxx::client conn{mongocxx::uri{}};

    double time_interval = 30;  //30 seconds
    int sliding_window_len = 6;
    int sliding_window_step = 1;
    // gen_datasets(conn, "sgtaxi_ref", "../density_test_dataset/sgtaxi.ds", sliding_window_len, sliding_window_step, time_interval);
    // gen_datasets(conn, "sgtaxi_afternoon_q", "../density_test_dataset/sgtaxi.q", sliding_window_len, sliding_window_step, time_interval, 0.01);

    time_interval = 600;  //30 seconds
    gen_datasets(conn, "interp_may_jun", "../density_test_dataset/vessel.ds", sliding_window_len, sliding_window_step, time_interval);
    gen_datasets(conn, "interp_aug_q", "../density_test_dataset/vessel.q", sliding_window_len, sliding_window_step, time_interval, 0.001);
}