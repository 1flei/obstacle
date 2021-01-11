#include <iostream>


#include "../util.h"
#include "interpolation.h"
#include "../mongodb/driver.hpp"

using bsoncxx::builder::basic::kvp;
using bsoncxx::builder::basic::make_array;
using bsoncxx::builder::basic::make_document;
// using bsoncxx::builder::stream::document;
// using bsoncxx::builder::stream::open_array;
// using bsoncxx::builder::stream::close_array;
// using bsoncxx::builder::stream::finalize;
using bsoncxx::to_json;

using namespace std;

const int sliding_window_len = 6;
const int sliding_window_step = 1;
const double time_interval = 600;
const double time_interval_freq = 30;


void test_interp_sorted()
{
    mongocxx::instance inst{};
    mongocxx::client conn{mongocxx::uri{}};

    cout.precision(16);

    auto collection = conn["traj_db"]["dsta_collection"];
    auto output_collection = conn["traj_db"]["interp_full"];
    {

        //first get all vessel mmsi
        auto mmsi_cursor = collection.distinct("tid", {});

        std::vector<std::string> all_mmsi;

        for(auto&& doc:mmsi_cursor){
            for(auto&& mmsi_utf8:doc["values"].get_array().value){
                std::string mmsi = mmsi_utf8.get_utf8().value.to_string();

                all_mmsi.push_back(std::move(mmsi));
            }
            // std::cout << to_json(doc) << std::endl;
        }

        fmt::print("all_mmsi.size()={}\n", all_mmsi.size());

        mongocxx::options::find opts;
        opts.sort(make_document(kvp("timestamp", 1)));
        // opts.limit(nRecordsToTest);
        for(const std::string& mmsi:all_mmsi){       
            int tid = get_id_from_mmsi(mmsi, conn);
            fmt::print("mmsi={}, tid={}\n", mmsi, tid);         
            auto maybe_doc = output_collection.find_one(make_document(kvp("tid", tid)));
            if(maybe_doc){
                //already exists
                fmt::print("already exists\n");
                continue;
            }

            auto cursor = collection.find(make_document(kvp("tid", mmsi)), opts);
            std::vector<Pnt3> cur_traj_vec;
            for(auto&& doc:cursor){
                Pnt3 p = doc_to_pnt3(doc);
                cur_traj_vec.emplace_back(std::move(p));
            }
           fmt::print("cur_traj_vec.size()={}\n", cur_traj_vec.size());
            //push
            auto interpolated_data = interpolate(cur_traj_vec, time_interval);
            if(!interpolated_data.empty()){
                auto ps = bsoncxx::builder::basic::array{};
                for(auto&& pp:interpolated_data){
                    ps.append(pp.t);
                    ps.append(pp.x);
                    ps.append(pp.y);
                }
                output_collection.insert_one(make_document(kvp("tid", tid), kvp("pnts", ps)));
            }
        }
    }
}

void test_interp_between(
    mongocxx::client& conn, 
    mongocxx::v_noabi::collection& collection, 
    mongocxx::v_noabi::collection& output_collection, 
    double t_start, 
    double t_end)
{
    cout.precision(16);
    {

        //first get all vessel mmsi
        auto mmsi_cursor = collection.distinct("tid", {});

        std::vector<std::string> all_mmsi;

        for(auto&& doc:mmsi_cursor){
            for(auto&& mmsi_utf8:doc["values"].get_array().value){
                std::string mmsi = mmsi_utf8.get_utf8().value.to_string();

                all_mmsi.push_back(std::move(mmsi));
            }
            // std::cout << to_json(doc) << std::endl;
        }

        fmt::print("all_mmsi.size()={}\n", all_mmsi.size());

        mongocxx::options::find opts;
        // opts.sort(make_document(kvp("timestamp", 1)));
        // opts.limit(nRecordsToTest);
        for(const std::string& mmsi:all_mmsi){       
            int tid = get_id_from_mmsi(mmsi, conn);
            fmt::print("mmsi={}, tid={}\n", mmsi, tid);         
            auto maybe_doc = output_collection.find_one(make_document(kvp("tid", tid)));
            if(maybe_doc){
                //already exists
                fmt::print("already exists\n");
                continue;
            }

            auto query_json = fmt::format("{{\"tid\": \"{}\", \"timestamp\": {{\"$gt\" : {}, \"$lt\": {} }}}}", mmsi, t_start, t_end);
            // fmt::print("{}\n", query_json);
            auto cursor = collection.find(bsoncxx::from_json(query_json), opts);
            std::vector<Pnt3> cur_traj_vec;
            for(auto&& doc:cursor){
                Pnt3 p = doc_to_pnt3(doc);
                cur_traj_vec.emplace_back(std::move(p));
            }
            fmt::print("found, sorting\n");
            std::sort(cur_traj_vec.begin(), cur_traj_vec.end(), [](const Pnt3& p0, const Pnt3& p1){
                return p0.t < p1.t;
            });
            fmt::print("cur_traj_vec.size()={}\n", cur_traj_vec.size());
            //push
            auto interpolated_data = interpolate(cur_traj_vec, time_interval);
            if(!interpolated_data.empty()){
                auto ps = bsoncxx::builder::basic::array{};
                for(auto&& pp:interpolated_data){
                    ps.append(pp.t);
                    ps.append(pp.x);
                    ps.append(pp.y);
                }
                output_collection.insert_one(make_document(kvp("tid", tid), kvp("pnts", ps)));
            }
        }
    }
}

//parse the traj_db.dsta_collection and interpolate and save to another collection

void test_interp(
    mongocxx::client& conn, 
    mongocxx::v_noabi::collection& collection, 
    mongocxx::v_noabi::collection& output_collection, 
    double t_start, 
    double t_end)
{

    cout.precision(16);
    std::unordered_map<std::string, std::vector<Pnt3> > tid_traj;

    {
        mongocxx::options::find opts;
        // opts.sort(make_document(kvp("tid", 1)));
        // opts.limit(nRecordsToTest);

        auto query_json = fmt::format("{{\"timestamp\": {{\"$gt\" : {}, \"$lt\": {} }}}}", t_start, t_end);
        auto cursor = collection.find(bsoncxx::from_json(query_json), opts);
        int cnt = 0;
        for(auto&& doc:cursor){
            auto[mmsi, p] = doc_to_pnt3_mmsi(doc);
            // std::cout << tid << ", " << lat << ", " << lng << ", " << timestamp << std::endl;
            tid_traj[mmsi].emplace_back(std::move(p));
            if(++cnt % 100000==0){
                fmt::print("cnt={}\n", cnt);
            }
        }
    }
    std::vector<bsoncxx::document::value> docs;
    for(auto&& p:tid_traj){
        std::sort(p.second.begin(), p.second.end(), [](const Pnt3& p0, const Pnt3& p1){
            return p0.t < p1.t;
        });
        auto interpolated_data = interpolate(p.second, time_interval);
        if(interpolated_data.empty()){
            continue;
        }
        auto ps = bsoncxx::builder::basic::array{};
        for(auto&& pp:interpolated_data){
            ps.append(pp.t);
            ps.append(pp.x);
            ps.append(pp.y);
        }
        int tid = get_id_from_mmsi(p.first, conn);
        docs.push_back(make_document(kvp("tid", tid), kvp("pnts", ps)));
    }
    output_collection.insert_many(docs);
}

void test_interp_bbox(
    mongocxx::client& conn, 
    mongocxx::v_noabi::collection& collection, 
    mongocxx::v_noabi::collection& output_collection, 
    const BBox& bb)
{

    cout.precision(16);
    std::unordered_map<std::string, std::vector<Pnt3> > tid_traj;

    {
        mongocxx::options::find opts;
        // opts.sort(make_document(kvp("tid", 1)));
        // opts.limit(nRecordsToTest);

        auto query_json = fmt::format("{{\"timestamp\": {{\"$gt\" : {}, \"$lt\": {} }}}}", bb.mint, bb.maxt);
        auto cursor = collection.find(bsoncxx::from_json(query_json), opts);
        int cnt = 0;
        for(auto&& doc:cursor){
            auto[mmsi, p] = doc_to_pnt3_mmsi(doc);
            if(++cnt % 100000==0){
                fmt::print("cnt={}\n", cnt);
            }
            if(p.x > bb.maxlng || p.x < bb.minlng || p.y > bb.maxlat || p.y < bb.minlat || p.t < bb.mint || p.t> bb.maxt) {
                continue;
            }
            // std::cout << tid << ", " << lat << ", " << lng << ", " << timestamp << std::endl;
            tid_traj[mmsi].emplace_back(std::move(p));
        }
    }
    std::vector<bsoncxx::document::value> docs;
    for(auto&& p:tid_traj){
        std::sort(p.second.begin(), p.second.end(), [](const Pnt3& p0, const Pnt3& p1){
            return p0.t < p1.t;
        });
        auto interpolated_data = interpolate(p.second, time_interval);
        if(interpolated_data.empty()){
            continue;
        }
        auto ps = bsoncxx::builder::basic::array{};
        for(auto&& pp:interpolated_data){
            ps.append(pp.t);
            ps.append(pp.x);
            ps.append(pp.y);
        }
        int tid = get_id_from_mmsi(p.first, conn);
        docs.push_back(make_document(kvp("tid", tid), kvp("pnts", ps)));
    }
    output_collection.insert_many(docs);
}


void test_interp_bbox_remove_anomaly(
    mongocxx::client& conn, 
    mongocxx::v_noabi::collection& collection, 
    mongocxx::v_noabi::collection& output_collection, 
    const BBox& bb, 
    double dist_constrint2)
{

    cout.precision(16);
    std::unordered_map<std::string, std::vector<Pnt3> > tid_traj;

    {
        mongocxx::options::find opts;
        // opts.sort(make_document(kvp("tid", 1)));
        // opts.limit(nRecordsToTest);

        auto query_json = fmt::format("{{\"timestamp\": {{\"$gt\" : {}, \"$lt\": {} }}}}", bb.mint, bb.maxt);
        auto cursor = collection.find(bsoncxx::from_json(query_json), opts);
        int cnt = 0;
        for(auto&& doc:cursor){
            auto[mmsi, p] = doc_to_pnt3_mmsi(doc);
            if(++cnt % 100000==0){
                fmt::print("cnt={}\n", cnt);
            }
            if(p.x > bb.maxlng || p.x < bb.minlng || p.y > bb.maxlat || p.y < bb.minlat || p.t < bb.mint || p.t> bb.maxt) {
                continue;
            }
            // std::cout << tid << ", " << lat << ", " << lng << ", " << timestamp << std::endl;
            tid_traj[mmsi].emplace_back(std::move(p));
        }
    }
    std::vector<bsoncxx::document::value> docs;
    for(auto&& p:tid_traj){
        std::sort(p.second.begin(), p.second.end(), [](const Pnt3& p0, const Pnt3& p1){
            return p0.t < p1.t;
        });
        auto interpolated_data = interpolate(p.second, time_interval, dist_constrint2);
        if(interpolated_data.empty()){
            continue;
        }
        auto ps = bsoncxx::builder::basic::array{};
        for(auto&& pp:interpolated_data){
            ps.append(pp.t);
            ps.append(pp.x);
            ps.append(pp.y);
        }
        int tid = get_id_from_mmsi(p.first, conn);
        docs.push_back(make_document(kvp("tid", tid), kvp("pnts", ps)));
    }
    output_collection.insert_many(docs);
}


void interp_daily(
    mongocxx::client& conn, 
    mongocxx::v_noabi::collection& collection, 
    mongocxx::v_noabi::collection& output_collection, 
    int hm_start, 
    int hm_end, 
    const std::string& given_state = "")
{
    cout.precision(16);
    std::unordered_map<std::string, std::vector<Pnt3> > tid_traj;

    {
        mongocxx::options::find opts;
        // opts.sort(make_document(kvp("tid", 1)));
        // opts.limit(nRecordsToTest);

        auto cursor = collection.find({}, opts);
        int cnt = 0;
        for(auto&& doc:cursor){
            auto[mmsi, state, p] = doc_to_pnt3_mmsi_state(doc);
            int timestamp = p.t;
            if(given_state!="" && given_state != state) {
                continue;
            }

            int hour_and_min = timestamp%(24*3600) / 60;

            // fmt::print("tid={}, p={}, p.t={}, timestamp={}, hour_and_min={}, hm_start={}, hm_end={}\n", mmsi, p, p.t, timestamp, hour_and_min, hm_start, hm_end);

            if(hour_and_min > hm_start && hour_and_min < hm_end){
                tid_traj[mmsi].emplace_back(std::move(p));
                if(++cnt % 100000==0){
                    fmt::print("cnt={}\n", cnt);
                }
            }
        }
    }
    std::vector<bsoncxx::document::value> docs;
    for(auto&& p:tid_traj){
        std::sort(p.second.begin(), p.second.end(), [](const Pnt3& p0, const Pnt3& p1){
            return p0.t < p1.t;
        });
        auto interpolated_data = interpolate(p.second, time_interval_freq);
        if(interpolated_data.empty()){
            continue;
        }
        auto ps = bsoncxx::builder::basic::array{};
        for(auto&& pp:interpolated_data){
            ps.append(pp.t);
            ps.append(pp.x);
            ps.append(pp.y);
        }
        int tid = get_id_from_mmsi(p.first, conn);
        docs.push_back(make_document(kvp("tid", tid), kvp("pnts", ps)));
    }
    output_collection.insert_many(docs);
}

void interp_ship()
{
    mongocxx::instance inst{};
    mongocxx::client conn{mongocxx::uri{}};
    
    const std::string colNameRef = "dsta_collection";
    const std::string colNameAug = "dsta_collection_aug";

    // auto collection = conn["traj_db"][collectionName];
    auto collectionRef = conn["traj_db"]["dsta_collection"];
    auto collectionAug = conn["traj_db"]["dsta_collection_aug"];

    // THORCO CLOUD (vessel 1)
    auto outputCollectionMay0 = conn["traj_db"]["interp_may0"];
    auto outputCollectionMay1 = conn["traj_db"]["interp_may1"];
    auto outputCollectionMayLarger = conn["traj_db"]["interp_may_larger"];
    auto outputCollectionAug = conn["traj_db"]["interp_aug"];

    auto outputCollectionMay_q = conn["traj_db"]["interp_may"];

    // Cai Jun 3 (vessel 2)
    auto outputCollectionJun_q = conn["traj_db"]["interp_jun_q"];

    //HARITA_BERLIAN
    auto outputCollectionAug_q = conn["traj_db"]["interp_aug_q"];

    double t_may = 1493568000;
    double t_jun = 1496246400;
    double t_jul = 1498838400;
    double t_aug = 1501516800;
    double t_sep = 1504195200;
    double dist_constraint2 = 0.05;
    // test_interp(conn, collectionRef, outputCollectionMay, t_may, t_jun);
    // may_q
    // test_interp_bbox(conn, collectionRef, outputCollectionMay0, BBox(103.60, 104.20, 1.0, 1.4, t_may, (t_may+t_jun)/2));
    
    // june_q
    // test_interp_bbox(conn, collectionRef, outputCollectionJun_q, BBox(104.1, 104.8, 1.25, 1.65, t_jun, (t_jun+t_jul)/2));
    // test_interp_bbox_remove_anomaly(conn, collectionRef, outputCollectionJun_q,  
    //     BBox(104.1, 104.8, 1.25, 1.65, t_jun, (t_jun+t_jul)/2), dist_constraint2);

    // aug_q 
    // test_interp_bbox(conn, collectionAug, outputCollectionAug_q, BBox(103.9, 104.5, 1.05, 1.45, t_aug, 0.75*t_aug+0.25*t_sep));
    // test_interp_bbox_remove_anomaly(conn, collectionAug, outputCollectionAug_q, 
    //     BBox(103.9, 104.5, 1.05, 1.45, t_aug, 0.75*t_aug+0.25*t_sep), dist_constraint2);


    // aug_ref
    // test_interp_between(conn, collectionAug, outputCollectionAug, t_aug, t_sep);
    test_interp_bbox_remove_anomaly(conn, collectionAug, outputCollectionAug, 
        BBox(-1e4, 1e4, -1e4, 1e4, t_aug, t_sep), dist_constraint2);

    // may_jun_ref
    // auto outputCollectionMayJunRef = conn["traj_db"]["interp_may_jun"];
    // // test_interp_bbox(conn, collectionRef, outputCollectionMayJunRef, BBox(-1e4, 1e4, -1e4, 1e4, t_may, t_jul) );
    // test_interp_bbox_remove_anomaly(conn, collectionRef, outputCollectionMayJunRef, 
    //     BBox(-1e4, 1e4, -1e4, 1e4, t_may, t_jul), dist_constraint2);
    // test_interp_between(conn, collectionRef, outputCollectionMayJunRef, t_may, t_jul);
}

void interp_sgtaxi()
{
    mongocxx::instance inst{};
    mongocxx::client conn{mongocxx::uri{}};
    
    const std::string colName = "sgtaxi";

    // auto collection = conn["traj_db"][collectionName];
    auto collection = conn["traj_db"]["sgtaxi"];

    auto outputCollectionRef = conn["traj_db"]["sgtaxi_ref"];

    auto outputCollectionMorningRef = conn["traj_db"]["sgtaxi_morning_ref"];
    auto outputCollectionMorningQ = conn["traj_db"]["sgtaxi_morning_q"];
    auto outputCollectionAfternoonRef = conn["traj_db"]["sgtaxi_afternoon_ref"];
    auto outputCollectionAfternoonQ = conn["traj_db"]["sgtaxi_afternoon_q"];

    int hm8 = 8*60;
    int hm10 = 10*60;
    int hm1930 = 19*60+30;
    int hm1730 = 17*60+30;

    // interp_daily(conn, collection, outputCollectionMorningRef, hm6am, hm8am, "FREE");
    interp_daily(conn, collection, outputCollectionRef, 0, 24*60, "FREE");
    // interp_daily(conn, collection, outputCollectionMorningQ, hm8am, hm9am, "FREE");
    interp_daily(conn, collection, outputCollectionMorningQ, hm8, hm10, "FREE");
    interp_daily(conn, collection, outputCollectionAfternoonQ, hm1730, hm1930, "FREE");
    // interp_daily(conn, collection, outputCollectionAfternoonQ, hm18pm, hm19pm, "FREE");
}


int main()
{
    interp_ship();
    // interp_sgtaxi();
    return 0;
}