#include "index.h"
#include "../mongodb/driver.hpp"
#include <iostream>
#include "my_pri_queue.h"

const int sliding_window_len = 6;
const int sliding_window_step = 1;
const double time_interval = 600.;

void build_from_database(mongocxx::client& conn, 
    hnswlib::HierarchicalNSW<double, SubtIdx>& appr_alg, 
    hnswlib::BruteforceSearch<double, SubtIdx>& bf_alg)
{
    auto ret = read_flat_trajs(conn);
    int cnt = 0;

    std::vector<Pnt3*> queries;

    for(auto&& p:ret){
        auto& t = p.second;
        auto& id = p.first;
        fmt::print("traj-{}, lens={}\n", id, t.get_num_pnts());
        sliding_window(t, sliding_window_len, sliding_window_step, time_interval, [&](const TrajView& tv){
            int offset = std::distance(&t[0], tv.p);

            appr_alg.addPoint(tv.p, SubtIdx(id, offset));
            bf_alg.addPoint(tv.p, SubtIdx(id, offset));
            ++cnt;
        });
    }
    fmt::print("cnt={}\n", cnt);
}


void index_test()
{
    mongocxx::instance inst{};
    mongocxx::client conn{mongocxx::uri{}};

    // HNSW hnsw;
    // hnswlib::L2Space;
    int maxObjs = 1000000;
    size_t M = 16;
    size_t ef_construction = 200;

    Pnt3DTWSpace dtwspace(sliding_window_len);
    hnswlib::HierarchicalNSW<double, SubtIdx> appr_alg(&dtwspace, maxObjs, M, ef_construction);
    hnswlib::BruteforceSearch<double, SubtIdx> bf_alg(&dtwspace, maxObjs);
    
    //since hnswlib reports error through exception
    try {
        appr_alg.loadIndex("interp.idx", &dtwspace, maxObjs);
        bf_alg.loadIndex("interp_bf.idx", &dtwspace);
    } catch(const std::exception& e) {
        fmt::print("exception when reading index: {}\n", e.what());
        fmt::print("so build from databases\n");

        build_from_database(conn, appr_alg, bf_alg);

        appr_alg.saveIndex("interp.idx");
        bf_alg.saveIndex("interp_bf.idx");
    }

    int nqueries = 100;
    auto traj_dict = read_flat_trajs(conn);
    std::string queryCollectionName = "interp_q";
    auto queries_dict = read_flat_trajs(conn, queryCollectionName);
    std::vector<TrajView> queries;

    for(auto&& p:queries_dict){
        int id = p.first;
        auto& t = p.second;
        sliding_window_n(t, nqueries-queries.size(), 
            sliding_window_len, sliding_window_step, time_interval, [&](const TrajView& tv){
            queries.push_back(tv);
        });
        if(queries.size()>=nqueries){
            break;
        }
    }

    fmt::print("start query\n");
    for(int i=0;i<queries.size();i++){
        fmt::print("query-{}; trajs={}\n", i, queries[i]);
        auto res0 = appr_alg.searchKnn(queries[i].p, 10);
        // auto res1 = bf_alg.searchKnn(queries[i].p, 10);

        while(!res0.empty()){
            auto& p0 = res0.top();
            fmt::print("    p0={}\n", p0);
            SubtIdx idx0 = p0.second;
            // fmt::print("        traj={}\n",TrajView(traj_dict[idx0.tid], idx0.offset, sliding_window_len));
            res0.pop();
        }
        auto res = searchKNN(appr_alg, queries[i].p, 10);
        fmt::print("res={}\n", res);
    }
}

void dtw_test()
{
    int ntrajs = 100;
    mongocxx::client conn{mongocxx::uri{}};
    auto traj_dict = read_flat_trajs(conn);
    std::string queryCollectionName = "interp_q";
    auto queries_dict = read_flat_trajs(conn, "interp_q");

    std::vector<TrajView> xs, ys;

    auto get_trajview = [&](const auto& dict, std::vector<TrajView>& ps){
        for(auto&& p:dict){
            int id = p.first;
            auto& t = p.second;
            sliding_window_n(t, ntrajs-ps.size(), 
            sliding_window_len, sliding_window_step, time_interval, [&](const TrajView& tv){
                ps.push_back(tv);
            });
            if(ps.size()>=ntrajs){
                break;
            }
        }
    };
    get_trajview(traj_dict, xs);
    get_trajview(queries_dict, ys);

    Pnt3DTWNormalizedSpace dtwspace(sliding_window_len);
    for(int i=0;i<ntrajs;i++){
        fmt::print("x={}\n", xs[i].to_flat_vec());
        fmt::print("y={}\n", ys[i].to_flat_vec());
        fmt::print("----dtw={}\n", dtwspace.get_dist_func()(xs[i].p, ys[i].p, dtwspace.get_dist_func_param()));
    }
}

void index_test2()
{
    mongocxx::instance inst{};
    mongocxx::client conn{mongocxx::uri{}};

    int maxObjs = 6000000;
    size_t M = 16;
    size_t ef_construction = 200;

    Pnt3DTWSpace dtwspace(sliding_window_len);
    hnswlib::HierarchicalNSW<double, SubtIdx> hnsw(&dtwspace, maxObjs, M, ef_construction);
    hnswlib::HierarchicalNSW<double, SubtIdx> hnsw_q(&dtwspace, maxObjs, M, ef_construction);


    std::string dataCollectionName = "interp_q";
    std::string queryCollectionName = "interp";
    std::string dataIndexFIle = fmt::format("{}.idx", dataCollectionName);
    std::string queryIndexFile = fmt::format("{}.idx", queryCollectionName);

    auto ref_dict = read_flat_trajs(conn, dataCollectionName);
    auto queries_dict = read_flat_trajs(conn, queryCollectionName);
    fmt::print("ref_dict.size()={}, queries_dict.size()={}\n", ref_dict.size(), queries_dict.size());   

    try {
        hnsw.loadIndex(dataIndexFIle, &dtwspace, maxObjs);
    } catch(const std::exception& e) {
        fmt::print("exception when reading index: {}\n", e.what());
        fmt::print("so build from databases\n");

        build_index(ref_dict, sliding_window_len, sliding_window_step, time_interval, hnsw);
        hnsw.saveIndex(dataIndexFIle);
    }

    try {
        //for test purpose
        hnsw_q.loadIndex(queryIndexFile, &dtwspace, maxObjs);
    } catch(const std::exception& e) {
        fmt::print("exception when reading index: {}\n", e.what());
        fmt::print("so build from databases\n");

        build_index(queries_dict, sliding_window_len, sliding_window_step, time_interval, hnsw_q);
        hnsw_q.saveIndex(queryIndexFile);
    }
    //build query index
    // build_index(queries_dict, hnsw_q);

    int numNeighbors = 10;
    int querycnt = 0;
    sliding_window_trajdict(queries_dict, sliding_window_len, sliding_window_step, time_interval, [&](const SubtIdx& qidx, const TrajView& tv){
        auto qvs = tv.to_flat_vec();
        fmt::print("----------------------------------\nquery-{}, qidx={}, trajs={}\n", querycnt, qidx, qvs);
        auto nr = searchKNN(hnsw, tv.p, numNeighbors);

        fmt::print("nr={}\n", nr);
        for(const auto& nr_pair:nr){
            double distq = nr_pair.first;
            SubtIdx idx = nr_pair.second;

            TrajView tv_refi(ref_dict[idx.tid], idx.offset, sliding_window_len);
            auto vs= tv_refi.to_flat_vec();

            double dtw_ = dtwspace.get_dist_func()(tv.p, tv_refi.p, dtwspace.get_dist_func_param());

            fmt::print("distance={}, dtw~={}; refi={}\n", distq, dtw_, vs);
        }
        if(++querycnt>100){
            exit(0);
        }
    });
}

int main()
{
    // index_test();
    dtw_test();
    // index_test2();
    return 0;
}