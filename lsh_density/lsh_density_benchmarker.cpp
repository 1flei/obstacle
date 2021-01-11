#include "boost/program_options.hpp"
#include <iostream>
#include "../my_timer.h"
#include "index.h"
#include "json.hpp"
#include "dataset_io.h"
#include "../lsh_density/lsh_density.h"
#include "../lsh_density/lsh_density_linear.h"

using namespace boost::program_options;


std::vector<double> brute_force_kde(const std::vector<std::vector<double> >& datasets, const std::vector<std::vector<double>>& queries, double sigma)
{
    std::vector<double> kdes;
    kdes.reserve(queries.size());

    for(auto& vec_q:queries) {
        double kde = 0;
        for(auto vec_di: datasets) {
            double dist2 = calc_l2_sqr<double>(vec_q.size(), &vec_q[0], &vec_di[0]);
            kde += exp(-dist2 * 0.5 / sigma/sigma);
        }
        kdes.push_back(kde);
    }
    return kdes;
}

void output_ground_truth_to_json(
    const std::vector<double>& ground_truth, 
    const std::vector<double>& estimated_density, 
    nlohmann::json& j_ret
)
{
    double mse = 0;
    double mre = 0.;
    for(int i=0;i<ground_truth.size();i++){
        fmt::print("q{}: gt={}, est_d={}\n", i, ground_truth[i], estimated_density[i]);
        mse += (ground_truth[i]-estimated_density[i]) *(ground_truth[i]-estimated_density[i]);

        mre += abs(ground_truth[i]-estimated_density[i]) /ground_truth[i];
    }
    j_ret["mse"] = mse;
    j_ret["mre"] = mre/ground_truth.size();
}

//parameters will be n_repeat, n_cm_repeat, ht_size
nlohmann::json test_lsh_denstiy_pie_cm(
    const std::vector<std::vector<double> >& datasets, 
    const std::vector<std::vector<double> >& queries, 
    const std::vector<double>& ground_truth, 
    const std::string& index_filename, 
    bool is_forcing_rebuilding_index, 
    int n_repeat, 
    int n_cm_repeat, 
    int ht_size, 
    double sigma,
    int k=8, int l=3, double r=5.12
)
{
    nlohmann::json j_ret;

    int dim = datasets[0].size();

    LSH_KDE index(dim, k, l, r*sigma, n_repeat, ht_size, n_cm_repeat);
    
    auto indexing_time = MyTimer::funcTime([&](){
        if(is_forcing_rebuilding_index){
            fmt::print("force building index from dataset\n");

            for(int i=0;i<datasets.size();i++){
                index.insert(datasets[i]);
            }
            if(index_filename!=""){
                index.saveIndex(index_filename);
            }
        } else{
            try {
                index.loadIndex(index_filename);
            } catch(const std::exception& e) {
                fmt::print("exception when reading index: {}\n", e.what());
                fmt::print("so build from databases\n");

                for(int i=0;i<datasets.size();i++){
                    index.insert(datasets[i]);
                }
                if(index_filename!=""){
                    index.saveIndex(index_filename);
                }
            }
        }
    });

    std::vector<double> estimated_density(queries.size());
    //query phase
    auto tot_query_time = MyTimer::funcTime([&](){
        for(int i=0;i<queries.size();i++) {
            estimated_density[i] = index.get_estimated_density(queries[i]);
        }
    });

    j_ret["indexing_time"] = indexing_time.count();
    j_ret["total_query_time"] = tot_query_time.count();
    j_ret["avg_query_time"] = tot_query_time.count() / queries.size();
    j_ret["n_repeat"] = n_repeat;
    j_ret["n_cm_repeat"] = n_cm_repeat;
    j_ret["ht_size"] = ht_size;
    // j_ret["est_densities"] = estimated_density;
    // j_ret["gt_densities"] = ground_truth;

    output_ground_truth_to_json(ground_truth, estimated_density, j_ret);

    j_ret["memory_usage"] = index.get_memory_usage();
    return j_ret;
}

//parameters will be n_repeat, n_cm_repeat, ht_size
nlohmann::json test_lsh_denstiy_linear(
    const std::vector<std::vector<double> >& datasets, 
    const std::vector<std::vector<double> >& queries, 
    const std::vector<double>& ground_truth, 
    const std::string& index_filename, 
    bool is_forcing_rebuilding_index, 
    int n_repeat, 
    int n_cm_repeat, 
    int ht_size, 
    double sigma
)
{
    nlohmann::json j_ret;

    int dim = datasets[0].size();

    LSH_KDE_LINEAR index(dim, 4*sigma, n_repeat, ht_size, n_cm_repeat);
    
    auto indexing_time = MyTimer::funcTime([&](){
        if(is_forcing_rebuilding_index){
            fmt::print("force building index from dataset\n");

            for(int i=0;i<datasets.size();i++){
                index.insert(datasets[i]);
            }
            if(index_filename!=""){
                index.saveIndex(index_filename);
            }
        } else{
            try {
                index.loadIndex(index_filename);
            } catch(const std::exception& e) {
                fmt::print("exception when reading index: {}\n", e.what());
                fmt::print("so build from databases\n");

                for(int i=0;i<datasets.size();i++){
                    index.insert(datasets[i]);
                }
                if(index_filename!=""){
                    index.saveIndex(index_filename);
                }
            }
        }
    });

    std::vector<double> estimated_density(queries.size());
    //query phase
    auto tot_query_time = MyTimer::funcTime([&](){
        for(int i=0;i<queries.size();i++) {
            estimated_density[i] = index.get_estimated_density(queries[i]);
        }
    });

    j_ret["indexing_time"] = indexing_time.count();
    j_ret["total_query_time"] = tot_query_time.count();
    j_ret["avg_query_time"] = tot_query_time.count() / queries.size();
    j_ret["n_repeat"] = n_repeat;
    j_ret["n_cm_repeat"] = n_cm_repeat;
    j_ret["ht_size"] = ht_size;
    // j_ret["est_densities"] = estimated_density;
    // j_ret["gt_densities"] = ground_truth;

    output_ground_truth_to_json(ground_truth, estimated_density, j_ret);
    j_ret["memory_usage"] = index.get_memory_usage();
    return j_ret;
}


int main(int argc, char **argv)
{
    srand(666);
    options_description desc("Allowed options");

    std::string dataset_filename, query_filename, ground_truth_filename, output_filename;

    double sigma = 0.001;
    int nqueries = 100;

    int ht_size = 4194301;
    int n_repeat = 64;
    int n_cm_repeat = 1;
    bool is_force_build_index = false;
    std::string index_filename;

    int bks_size;

    desc.add_options()
        ("help,h", "help message")
        ("dataset_filename,D", value(&dataset_filename), "data file name")
        ("query_filename,Q", value(&query_filename), "query file name")
        ("ground_truth_filename,G", value(&ground_truth_filename), "ground_truth_filename")
        ("output_filename,O", value(&output_filename), "output file name (usually xxx.json)")

        ("nqueries,q", value(&nqueries), "# queries")

        ("compute_ground_truth", "compute ground truth if specified")
        ("sigma", value(&sigma), "sigma for KDE")


        ("index_filename,I", value(&index_filename), "index file name")
        ("ht_size", value(&ht_size), "size of hash table for pie+cm")
        ("n_repeat", value(&n_repeat), "# tiems repeated")
        ("n_cm_repeat", value(&n_cm_repeat), "# times repeated for count min sketch")

        ("bks_size", value(&bks_size), "size of bottom-k sketch")

        ("lsh_density_pie", "using principle of inclusion and exclusion + count min sketch")
        ("lsh_density_linear", "using linaer regression")

        ("force_build_index", "force rebuild index")
    ;
    variables_map vm;
    store(parse_command_line(argc, argv, desc), vm);
    notify(vm);  

    if (vm.count("help")) {
        std::cout << desc << std::endl;
        return 1;
    }

    if (vm.count("force_build_index")) {
        is_force_build_index = true;
    }

    std::vector<std::vector<double > > datasets = read_dataset(dataset_filename);
    std::vector<std::vector<double > > queries = read_dataset(query_filename, nqueries);

    fmt::print("datasets.size()={}, queries.size()={}\n", datasets.size(), queries.size());

    //use lambda to apply RVO
    std::vector<double> gt = [&](){
        if(vm.count("compute_ground_truth")) {
            MyTimer::pusht();
            auto gt = brute_force_kde(datasets, queries, sigma);
            double time = MyTimer::popt();

            fmt::print("brute_force_time={}\n", time);
            write_ground_truth(ground_truth_filename, gt);
            return gt;
        } else{
            return read_ground_truth(ground_truth_filename, nqueries);
        }
    }();

    fmt::print("ground_truth.size()={}\n", gt.size());

    std::ofstream outf(output_filename, std::ofstream::app);

    // for(auto p:vm){
    //     fmt::print("{}\n", p.first);
    // }
    
    if(vm.count("lsh_density_pie")) {
        fmt::print("lsh_density_pie!!\n");
        auto res_json = test_lsh_denstiy_pie_cm(datasets, queries, gt, index_filename, is_force_build_index, n_repeat, n_cm_repeat, ht_size, sigma);
        outf << res_json.dump(4) << std::endl;
    }

    if(vm.count("lsh_density_linear")) {
        fmt::print("lsh_density_linear!!\n");
        auto res_json = test_lsh_denstiy_linear(datasets, queries, gt, index_filename, is_force_build_index, n_repeat, n_cm_repeat, ht_size, sigma);
        outf << res_json.dump(4) << std::endl;
    }
    return 0;
}