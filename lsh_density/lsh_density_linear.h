#pragma once

#include "e2eigen.h"
#include "tabu_hash.h"
#include <unordered_map>
#include <vector>
#include <fstream>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/unordered_map.hpp>

#ifdef _MSC_VER
#  include <intrin.h>
#  define __builtin_popcount __popcnt
#endif

//using inclusion-exclusion principle using flat hash
class LSH_KDE_LINEAR
{
public:
    int dim, nrepeat, ht_size, count_min_sketch_nrepeat;
    double r;
    std::vector<E2Eigen> lsh_hashers;
    // std::vector<std::unordered_map<E2Eigen::SigType, int32_t>> hash_tables;
    std::vector<std::vector<int32_t>> hash_tables;
    // will have (2**l -1)*nrepeat hash tables
    // each has size ht_size (which would better be a prime number)
    // utimately there will be a best combination for nrepeat and ht_size so that nrepeat*ht_size = constant

    std::vector<TabuHasher64> pert_hashers;
    //use tabulation hash to combine different hash values and generate count_min_sketch_nrepeat different hash values  
    //out of each lsh value
    //since each time it might need combine maximum l hashers
    //so the total size will be l*count_min_sketch_nrepeat

    // static inline std::array<double, 7> weights = {-0.09107084,  0.41073369,  0.41073369,  0.41073369,  0.41073369,
    //         -0.28826238, -0.28826238};
    // static inline std::array<int, 7> ks = {1, 4, 4, 4, 4, 15, 15};
    // static inline int ks_sum = 47;
    static inline std::array<double, 7> weights = {0.03581134, -0.87357163,  1.97736006,  1.97736006, -1.65327229,
       -1.6066025 ,  1.1412126};
    static inline std::array<int, 7> ks = {1, 3, 5, 5, 10, 11, 15};
    static inline int ks_sum = 50;
    static inline int nterms = 7;
public:

    LSH_KDE_LINEAR(int dim, double r, int nrepeat = 1, int ht_size = 4194301, int count_min_sketch_nrepeat = 1)
        : dim(dim), nrepeat(nrepeat), ht_size(ht_size), count_min_sketch_nrepeat(count_min_sketch_nrepeat), hash_tables(nterms*nrepeat)
    {
        lsh_hashers.reserve(nrepeat);
        for(int i=0;i<nrepeat;i++){
            lsh_hashers.emplace_back(dim, ks_sum, r, 1, i*E2Eigen::DEFAULT_SEED);
        }

        for(int i=0;i<hash_tables.size();i++){
            hash_tables[i].resize(ht_size);
        }

        pert_hashers.reserve(count_min_sketch_nrepeat*ks_sum);
        for(int i=0;i<count_min_sketch_nrepeat*ks_sum;i++){
            pert_hashers.emplace_back(i*TabuHasher32::DEFAULT_SEED);
        }
    }

    void insert(const std::vector<double> &x)
    {
        // fmt::print("x={}\n", x);
        for(int i=0; i<nrepeat; i++) {
            auto& lsh_hasher = lsh_hashers[i];
            //sig = 47 hash values
            auto sig = lsh_hasher.get_sig_concatenated(&x[0]);

            for(int k=0;k<count_min_sketch_nrepeat;k++){
                int cur_hasher_idx = 0;
                for(int j=0;j<nterms;j++){
                    int kj = ks[j];
                    uint64_t hash_value = 0;
                    for(int jj=0;jj<kj;jj++){
                        hash_value = hash_value ^ pert_hashers[cur_hasher_idx].hash(sig[cur_hasher_idx]);
                        ++cur_hasher_idx;
                    }
                    int hash_table_idx = i*nterms + j;
                    hash_tables[hash_table_idx][hash_value%ht_size] += 1;
                }
                // assert(cur_hasher_idx==ks_sum);
            }
        }
    }

    double get_estimated_density(const std::vector<double> &x)
    {
        // fmt::print("x={}\n", x);
        std::vector<double> cnts(nrepeat*(nterms)*count_min_sketch_nrepeat);

        for(int i=0; i<nrepeat; i++) {
            auto& lsh_hasher = lsh_hashers[i];
            auto sig = lsh_hasher.get_sig_concatenated(&x[0]);
            // fmt::print("for hasher-{}, sig={}\n", i, sig);
            for(int k=0;k<count_min_sketch_nrepeat;k++){
                int cur_hasher_idx = 0;
                for(int j=0;j<nterms;j++){
                    int kj = ks[j];
                    uint64_t hash_value = 0;
                    for(int jj=0;jj<kj;jj++){
                        hash_value = hash_value ^ pert_hashers[cur_hasher_idx].hash(sig[cur_hasher_idx]);
                        ++cur_hasher_idx;
                    }

                    int hash_table_idx = i*nterms + j;
                    int cnt_idx = (i*(nterms)+j)*count_min_sketch_nrepeat+k;
                    // hash_tables[hash_table_idx][hash_value%ht_size] += 1;
                    cnts[cnt_idx] = hash_tables[hash_table_idx][hash_value%ht_size];
                }
                // assert(cur_hasher_idx==ks_sum);
            }
        }

        //mean over nrepeat
        double sum = 0;
        for(int i=0;i<nrepeat;i++){
            //sum the cnt up using the exclusion and inclusion principle
            double current_density_estimator = 0;
            for(int j=0;j<nterms;j++) {
                //count-min sketch like estimator
                int mincnt = 1e8;
                for(int k=0;k<count_min_sketch_nrepeat;k++){
                    int idx = (i*(nterms)+j)*count_min_sketch_nrepeat+k;
                    if(cnts[idx] < mincnt){
                        mincnt = cnts[idx];
                    }
                    // printf("i=%d, j=%d, k=%d, cnt=%f\n", i, j, k, cnts[idx]);
                }
                current_density_estimator += weights[j] * mincnt;
            }
            sum += current_density_estimator;
        }
        return sum/nrepeat;

        // std::nth_element(cnts.begin(), cnts.begin() + nrepeat / 2, cnts.end()); 
        // return cnts[nrepeat / 2];
        //return average
        // double sum = 0;
        // for(int i=0;i<nrepeat;i++){
        //     sum += cnts[i];
        // }
        // return sum/ nrepeat;
    }

    template <class Archive>
    void serialize(Archive &ar, const unsigned int )
    {
        ar & dim;
        ar & nrepeat;
        ar & ht_size;
        ar & count_min_sketch_nrepeat;
        ar & lsh_hashers;
        ar & hash_tables;
        ar & pert_hashers;
    }

    void saveIndex(const std::string& filename)
    {
        std::ofstream f(filename);
        boost::archive::binary_oarchive oa(f);
        oa & (*this);
    }

    void loadIndex(const std::string& filename)
    {
        std::ifstream f(filename);
        boost::archive::binary_iarchive oa(f);
        oa & (*this);
    }

    int64_t get_memory_usage()
    {
        int64_t ret = sizeof(*this);
        //std::vector<E2Eigen> lsh_hashers;
        for(int i=0;i<lsh_hashers.size();i++){
            ret += lsh_hashers[i].get_memory_usage();
        }
        // std::vector<std::vector<int32_t>> hash_tables;
        for(int i=0;i<hash_tables.size();i++){
            ret += sizeof(hash_tables[i]);
            ret += hash_tables[i].capacity()*sizeof(int32_t);
        }
        // std::vector<TabuHasher32> pert_hashers;
        for(int i=0;i<pert_hashers.size();i++){
            ret += pert_hashers[i].get_memory_usage();
        }
        return ret;
    }
};

