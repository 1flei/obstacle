#include <iostream>


#include "../util.h"
#include "interpolation.h"
#include "../mongodb/driver.hpp"
#include <unordered_map>
#include <unordered_set>
#include "json.hpp"
#include <fstream>
#include "../index/subtidx.h"

using bsoncxx::builder::basic::kvp;
using bsoncxx::builder::basic::make_array;
using bsoncxx::builder::basic::make_document;
using bsoncxx::to_json;

using namespace std;

//parameters for ship dataset
const int sliding_window_len = 6;
const int sliding_window_step = 1;
const double time_interval = 600;

//tata_cable heatmap
void tatacable_heatmap()
{
    std::unordered_map<std::tuple<int, int>, int > heatmap;
    mongocxx::instance inst{};
    mongocxx::client conn{mongocxx::uri{}};

    cout.precision(16);

    int cnt = 0;

    double grid_width = 0.0002; //0.29/0.0002 = 1450

    //filter by me
    auto collection = conn["traj_db"]["dsta_collection"];
    {
        mongocxx::options::find opts;
        // opts.sort(make_document(kvp("tid", 1)));
        // opts.limit(nRecordsToTest);

        auto cursor = collection.find({}, opts);
        for(auto&& doc:cursor){
            // [1.0, 1.2] \times [103.42, 103.71]
            double lat = doc["lat"].get_double();
            double lng = doc["lng"].get_double();
            cnt++;
            if(cnt%100000==0){
                std::cout << cnt << std::endl;
            }
            if(lng > 103.71 || lng < 103.42 || lat > 1.2 || lat < 1.0) {
                continue;
            }
            int lat_grid = floor(lat/grid_width+0.5);
            int lng_grid = floor(lng/grid_width+0.5);

            heatmap[make_tuple(lng_grid, lat_grid)] += 1;
            // std::cout << to_json(doc) << std::endl;
        }
    }

    using nlohmann::json;
    json j = heatmap;
    std::cout << j.dump() << std::endl;
}
//tata_cable heatmap
void tatacable_heatmap_aug()
{
    std::unordered_map<std::tuple<int, int>, int > heatmap;
    mongocxx::instance inst{};
    mongocxx::client conn{mongocxx::uri{}};

    cout.precision(16);

    int cnt = 0;

    double grid_width = 0.0002; //0.29/0.0002 = 1450

    //filter by me
    auto collection = conn["traj_db"]["dsta_collection_aug"];
    {
        mongocxx::options::find opts;
        // opts.sort(make_document(kvp("tid", 1)));
        // opts.limit(nRecordsToTest);

        auto cursor = collection.find({}, opts);
        for(auto&& doc:cursor){
            // [1.0, 1.2] \times [103.42, 103.71]
            double lat = doc["lat"].get_double();
            double lng = doc["lng"].get_double();
            double timestamp = doc["timestamp"].get_double();
            cnt++;
            if(cnt%100000==0){
                std::cout << cnt << std::endl;
            }
            if(lng > 103.71 || lng < 103.42 || lat > 1.2 || lat < 1.0) {
                continue;
            }
            //Aug 1st - Aug 15th
            if(timestamp < 1501516801 || timestamp> 1502812800) {
                continue;
            }
            int lat_grid = floor(lat/grid_width+0.5);
            int lng_grid = floor(lng/grid_width+0.5);

            heatmap[make_tuple(lng_grid, lat_grid)] += 1;
            // std::cout << to_json(doc) << std::endl;
        }
    }

    using nlohmann::json;
    json j = heatmap;
    std::cout << j.dump() << std::endl;
}


void sunken_vessel_HARITA_BERLIAN()
{
    // [104.10, 104.25] * [1.25, 1.35]
    std::unordered_map<std::tuple<int, int>, int > heatmap;
    mongocxx::instance inst{};
    mongocxx::client conn{mongocxx::uri{}};

    cout.precision(16);

    int cnt = 0;

    double grid_width = 0.0002; //0.29/0.0002 = 1450

    //filter by me
    auto collection = conn["traj_db"]["dsta_collection"];
    {
        mongocxx::options::find opts;
        // opts.sort(make_document(kvp("tid", 1)));
        // opts.limit(nRecordsToTest);

        auto cursor = collection.find({}, opts);
        for(auto&& doc:cursor){
            // [1.0, 1.2] \times [103.42, 103.71]
            double lat = doc["lat"].get_double();
            double lng = doc["lng"].get_double();
            cnt++;
            if(cnt%100000==0){
                std::cout << cnt << std::endl;
            }
            if(lng > 104.25 || lng < 104.10 || lat > 1.35 || lat < 1.25) {
                continue;
            }
            int lat_grid = floor(lat/grid_width+0.5);
            int lng_grid = floor(lng/grid_width+0.5);

            heatmap[make_tuple(lng_grid, lat_grid)] += 1;
            // std::cout << to_json(doc) << std::endl;
        }
    }

    using nlohmann::json;
    json j = heatmap;
    std::cout << j.dump() << std::endl;
}


void sunken_vessel_HARITA_BERLIAN_aug()
{
    // [104.10, 104.25] * [1.25, 1.35]
    std::unordered_map<std::tuple<int, int>, int > heatmap;
    mongocxx::instance inst{};
    mongocxx::client conn{mongocxx::uri{}};

    cout.precision(16);

    int cnt = 0;

    double grid_width = 0.0002; //0.29/0.0002 = 1450

    //filter by me
    auto collection = conn["traj_db"]["dsta_collection_aug"];
    {
        mongocxx::options::find opts;
        // opts.sort(make_document(kvp("tid", 1)));
        // opts.limit(nRecordsToTest);

        auto cursor = collection.find({}, opts);
        for(auto&& doc:cursor){
            // [1.0, 1.2] \times [103.42, 103.71]
            double lat = doc["lat"].get_double();
            double lng = doc["lng"].get_double();
            double timestamp = doc["timestamp"].get_double();
            cnt++;
            if(cnt%100000==0){
                std::cout << cnt << std::endl;
            }
            if(lng > 104.25 || lng < 104.10 || lat > 1.35 || lat < 1.25) {
                continue;
            }
            //Aug 1st to 15th
            if(timestamp < 1501516801 || timestamp> 1502812800) {
                continue;
            }
            int lat_grid = floor(lat/grid_width+0.5);
            int lng_grid = floor(lng/grid_width+0.5);

            heatmap[make_tuple(lng_grid, lat_grid)] += 1;
            // std::cout << to_json(doc) << std::endl;
        }
    }

    using nlohmann::json;
    json j = heatmap;
    std::cout << j.dump() << std::endl;
}


void get_heatmap_json(
        mongocxx::v_noabi::collection& collection, 
        const BBox& bb, 
        const std::string& out_filename, 
        double grid_width = 0.0002)
{
    fmt::print("{}\n", out_filename);
    // [103. 80, 104. 00] * [1.15, 1.25] 
    std::unordered_map<std::tuple<int, int>, int > heatmap;

    cout.precision(16);

    int cnt = 0;
    int hcnt = 0;

    {
        mongocxx::options::find opts;
        // opts.sort(make_document(kvp("tid", 1)));
        // opts.limit(nRecordsToTest);

        auto query_json = fmt::format("{{\"timestamp\": {{\"$gt\" : {}, \"$lt\": {} }}}}", bb.mint, bb.maxt);
        std::cout << query_json << std::endl;
        auto query_doc = bsoncxx::from_json(query_json);

        auto cursor = collection.find(query_doc.view(), opts);
        for(auto&& doc:cursor){
            // [1.0, 1.2] \times [103.42, 103.71]
            Pnt3 p = doc_to_pnt3(doc);
            double lng = p.x;
            double lat = p.y;
            double timestamp = p.t;
            // fmt::print(" {}, {}, {} | {}, {}, {}, {}, {}, {}\n", lng, lat, timestamp, bb.maxlng, bb.minlng, bb.maxlat, bb.minlat, bb.mint, bb.maxt);
            cnt++;
            if(cnt%100000==0){
                std::cout << cnt << std::endl;
            }
            if(p.x > bb.maxlng || p.x < bb.minlng || p.y > bb.maxlat || p.y < bb.minlat || timestamp < bb.mint || timestamp> bb.maxt) {
                continue;
            }
            //Aug 1st to 15th
            int lat_grid = floor(lat/grid_width+0.5);
            int lng_grid = floor(lng/grid_width+0.5);

            heatmap[make_tuple(lng_grid, lat_grid)] += 1;
            // std::cout << to_json(doc) << std::endl;
            hcnt++;
            if(hcnt%10000==0){
                std::cout << "~" << hcnt << std::endl;
            }
        }
    }

    nlohmann::json j = heatmap;
    ofstream outf(out_filename);
    outf << j.dump();
}

void heatmap_from_interped_data(
        mongocxx::client& conn, 
        const std::string& queryCollectionName, 
        const BBox& bb, 
        const std::string& out_filename, 
        double grid_width = 0.0002)
{
    fmt::print("{}\n", out_filename);
    // [103. 80, 104. 00] * [1.15, 1.25] 
    std::unordered_map<std::tuple<int, int>, int > heatmap;

    cout.precision(16);

    int cnt = 0;
    int hcnt = 0;
    auto queries_dict = read_flat_trajs(conn, queryCollectionName);

    sliding_window_trajdict(queries_dict, sliding_window_len, sliding_window_step, time_interval, [&](const SubtIdx& qidx, const TrajView& tv){
        for(Pnt3& p:tv){
            double lng = p.x;
            double lat = p.y;
            double timestamp = p.t;
            cnt++;
            if(cnt%100000==0){
                std::cout << cnt << std::endl;
            }
            if(p.x > bb.maxlng || p.x < bb.minlng || p.y > bb.maxlat || p.y < bb.minlat || timestamp < bb.mint || timestamp> bb.maxt) {
                continue;
            }
            //Aug 1st to 15th
            int lat_grid = floor(lat/grid_width+0.5);
            int lng_grid = floor(lng/grid_width+0.5);

            heatmap[make_tuple(lng_grid, lat_grid)] += 1;
            // std::cout << to_json(doc) << std::endl;
            hcnt++;
            if(hcnt%10000==0){
                std::cout << "~" << hcnt << std::endl;
            }
        }
    });

    nlohmann::json j = heatmap;
    ofstream outf(out_filename);
    outf << j.dump();
}

void outputRawTrajJson(
    mongocxx::client& conn, 
    const std::string& collectionName, 
    const BBox& bb, 
    const std::string& output_filename
)
{
    fmt::print("{}\n", output_filename);
    // [103. 80, 104. 00] * [1.15, 1.25] 
    std::unordered_map<std::tuple<int, int>, int > heatmap;

    cout.precision(16);

    int cnt = 0;
    int hcnt = 0;
    auto queries_dict = read_flat_trajs(conn, collectionName);

    nlohmann::json j;
    sliding_window_trajdict(queries_dict, sliding_window_len, sliding_window_step, time_interval, [&](const SubtIdx& qidx, const TrajView& tv){
        auto vs = tv.to_flat_vec();
        j.push_back(vs);
    });

    ofstream outf(output_filename);
    outf << j.dump();
}

void get_heatmap_json_daily(
        mongocxx::v_noabi::collection& collection, 
        const BBox& bb, 
        const std::string& out_filename, 
        double grid_width = 0.0002, 
        int hourmin_start = 0, 
        int hourmin_end = 24*60, 
        const std::string& given_state = "")
{
    fmt::print("{}\n", out_filename);
    // [103. 80, 104. 00] * [1.15, 1.25] 
    std::unordered_map<std::tuple<int, int>, int > heatmap;

    cout.precision(16);

    int cnt = 0;
    int hcnt = 0;

    {
        mongocxx::options::find opts;
        auto cursor = collection.find({}, opts);
        for(auto&& doc:cursor){
            // [1.0, 1.2] \times [103.42, 103.71]
            auto[mmsi, state, p] = doc_to_pnt3_mmsi_state(doc);
            int timestamp = p.t;
            if(given_state!="" && given_state != state) {
                continue;
            }

            int hour_and_min = timestamp%(24*3600) / 60;


            if(++cnt % 1000000==0){
                fmt::print("cnt={}\n", cnt);
            }
            if(hour_and_min < hourmin_start || hour_and_min > hourmin_end){
                continue;
            }
            double lng = p.x;
            double lat = p.y;

            if(p.x > bb.maxlng || p.x < bb.minlng || p.y > bb.maxlat || p.y < bb.minlat) {
                continue;
            }
            // fmt::print("tid={}, p={}, p.t={}, timestamp={}, hour_and_min={}, hm_start={}, hm_end={}\n", mmsi, p, p.t, timestamp, hour_and_min, hourmin_start, hourmin_end);

            //Aug 1st to 15th
            int lat_grid = floor(lat/grid_width+0.5);
            int lng_grid = floor(lng/grid_width+0.5);

            heatmap[make_tuple(lng_grid, lat_grid)] += 1;
            // std::cout << to_json(doc) << std::endl;
            hcnt++;
            if(hcnt%100000==0){
                std::cout << "~" << hcnt << std::endl;
            }
        }
    }

    nlohmann::json j = heatmap;
    ofstream outf(out_filename);
    outf << j.dump();
}

int main()
{
    // test_interp();
    // sunken_vessel_HARITA_BERLIAN();
    // sunken_vessel_HARITA_BERLIAN_aug();
    mongocxx::instance inst{};
    mongocxx::client conn{mongocxx::uri{}};
    
    const std::string colNameRef = "dsta_collection";
    const std::string colNameAug = "dsta_collection_aug";

    // auto collection = conn["traj_db"][collectionName];
    auto collectionRef = conn["traj_db"]["dsta_collection"];
    auto collectionAug = conn["traj_db"]["dsta_collection_aug"];

    double t_may = 1493568000;
    double t_jun = 1496246400;
    double t_jul = 1498838400;
    double t_aug = 1501516800;
    double t_sep = 1504195200;

    //vessel Harita Berlian
    get_heatmap_json(collectionRef, BBox(104.10, 104.25, 1.25, 1.35, t_may, t_jul), "sunken_vessel_0_ref.json");
    get_heatmap_json(collectionAug, BBox(104.10, 104.25, 1.25, 1.35, t_aug, (t_aug+t_sep)/2), "sunken_vessel_0_q.json");

    //vessel THORCO CLOUD
    // get_heatmap_json(collectionRef, BBox(103.80, 104.00, 1.15, 1.25, t_may, (t_may+t_jun)/2), "sunken_vessel_1_may.json");
    // get_heatmap_json(collectionAug, BBox(103.80, 104.00, 1.15, 1.25, t_aug, (t_aug+t_sep)/2), "sunken_vessel_1_aug.json");
    // heatmap_from_interped_data(conn, "interp_may1", BBox(103.80, 104.00, 1.15, 1.25, t_may, (t_may+t_jun)/2), "interp_may1_hm.json");
    // heatmap_from_interped_data(conn, "interp_aug", BBox(103.80, 104.00, 1.15, 1.25, t_aug, (t_aug+t_sep)/2), "interp_aug_hm.json");
    // outputRawTrajJson(conn, "interp_may1", BBox(103.80, 104.00, 1.15, 1.25, t_may, (t_may+t_jun)/2), "interp_may1_raw_traj.json");
    // outputRawTrajJson(conn, "interp_aug", BBox(103.80, 104.00, 1.15, 1.25, t_may, (t_may+t_jun)/2), "interp_aug_raw_traj.json");
    
    //vessel Cai Jun 3
    // get_heatmap_json(collectionRef, BBox(104.35, 104.55, 1.35, 1.5, t_jun, t_jul), "sunken_vessel_2_jun.json");
    // get_heatmap_json(collectionAug, BBox(104.35, 104.55, 1.35, 1.5, t_aug, t_sep), "sunken_vessel_2_aug.json");
    
    //vessel Putri Sea
    // get_heatmap_json(collectionRef, BBox(104.05, 104.25, 1.2, 1.35, t_may, t_jul), "sunken_vessel_3_may.json");
    // get_heatmap_json(collectionAug, BBox(104.05, 104.25, 1.2, 1.35, t_aug, t_sep), "sunken_vessel_3_aug.json");



    auto sgtaxi_collection_ref = conn["traj_db"]["sgtaxi"];
    double grid_width = 0.0002;
    int hm8 = 8*60;
    int hm10 = 10*60;
    int hm1930 = 19*60+30;
    int hm1730 = 17*60+30;
    //ERP
    // get_heatmap_json_daily(sgtaxi_collection_ref, BBox(103.6, 104, 1.23, 1.47, t_may, t_jul), "sgtaxi_heatmap.json", 
    //     grid_width, 0*60, 24*60, "FREE");
    // get_heatmap_json_daily(sgtaxi_collection_ref, BBox(103.6, 104, 1.23, 1.47, t_may, t_jul), "sgtaxi_heatmap_morning.json", 
    //     grid_width, hm8, hm10, "FREE");
    // get_heatmap_json_daily(sgtaxi_collection_ref, BBox(103.6, 104, 1.23, 1.47, t_may, t_jul), "sgtaxi_heatmap_afternoon.json", 
    //     grid_width, hm1730, hm1930, "FREE");

    return 0;
}