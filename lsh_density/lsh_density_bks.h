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

#include "bottom_k_sketch.h"

//try to implement a simple version of
//Fair Near Neighbor Search: Independent Range Sampling in High Dimensions

class LSH_KDE_BKS
{
public:
    using BKS = BottomKSketch32;
    using BKS_SET = BottomKSketch32::Set;

    int dim, k, l, nrepeat;
    std::vector<E2Eigen> lsh_hashers;
    std::vector<BKS> bks_managers;
    std::vector<std::unordered_map<E2Eigen::SigType, BKS_SET>> hash_tables;

public:

    LSH_KDE_BKS(int dim, int k, int l, double r, int bottom_k_sketch_size = 64, int nrepeat = 1)
        : dim(dim), k(k), l(l), nrepeat(nrepeat), hash_tables(l*nrepeat)
    {
        bks_managers.reserve(nrepeat);
        for (int i = 0; i < nrepeat; i++)
        {
            //ideally will be the seed generated from the seed of LSH_KDE_BKS
            //TODO: will add this feature later on
            bks_managers.emplace_back(bottom_k_sketch_size, i);
        }

        lsh_hashers.reserve(nrepeat);
        for(int i=0;i<nrepeat;i++){
            lsh_hashers.emplace_back(dim, k*l, r, k, i*E2Eigen::DEFAULT_SEED);
        }
    }

    void insert(const std::vector<double> &x, int id)
    {
        for(int i=0; i<nrepeat; i++) {
            auto& lsh_hasher = lsh_hashers[i];
            auto& bks_manager = bks_managers[i];
            auto sig = lsh_hasher.get_sig_concatenated(&x[0]);
            for(int j=0; j<l; j++){
                int hash_table_idx = i*l + j;
                auto& hash_table = hash_tables[hash_table_idx];
                
                bks_manager.insert(id, hash_table[sig[j]]);
            }
        }
    }

    double get_estimated_density(const std::vector<double> &x)
    {
        std::vector<double> cnts(nrepeat);
        for(int i=0;i<nrepeat;i++){
            auto& lsh_hasher = lsh_hashers[i];
            auto& bks_manager = bks_managers[i];
            auto sig = lsh_hasher.get_sig_concatenated(&x[0]);

            BKS_SET sigset;
            for(int j=0;j<l;j++){
                int hash_table_idx = i*l + j;
                auto& hash_table = hash_tables[hash_table_idx];
                bks_manager.merge(sigset, hash_table[sig[j]]);
            }
            cnts[i] = bks_manager.get_cardinality(sigset);
        }

        // std::nth_element(cnts.begin(), cnts.begin() + nrepeat / 2, cnts.end()); 
        // return cnts[nrepeat / 2];
        //return average
        double sum = 0;
        for(int i=0;i<nrepeat;i++){
            sum += cnts[i];
        }
        return sum/ nrepeat;
    }

    template <class Archive>
    void serialize(Archive &ar, const unsigned int )
    {
        ar & dim;
        ar & k;
        ar & l;
        ar & nrepeat;
        ar & lsh_hashers;
        ar & bks_managers;
        ar & hash_tables;
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
        // std::vector<E2Eigen> lsh_hashers;
        int64_t ret = sizeof(*this);
        //std::vector<E2Eigen> lsh_hashers;
        for(int i=0;i<lsh_hashers.size();i++){
            ret += lsh_hashers[i].get_memory_usage();
        }
        // std::vector<BKS> bks_managers;
        for(int i=0;i<bks_managers.size();i++){
            ret += bks_managers[i].get_memory_usage();
        }
        // std::vector<std::unordered_map<E2Eigen::SigType, BKS_SET>> hash_tables;
        for(int i=0;i<hash_tables.size();i++){
            ret += sizeof(hash_tables[i]);
            int64_t bucket_size = sizeof(void*);
            ret += (int64_t)hash_tables[i].bucket_count() * bucket_size;
            for(auto&& p:hash_tables[i]){
                ret += p.second.get_memory_usage() + sizeof(void*) + sizeof(E2Eigen::SigType);
                // fmt::print("3.5.5!!!!!ret={}, step={}\n", ret, p.second.get_memory_usage() + sizeof(void*) + sizeof(E2Eigen::SigType));
            }
        }
        return ret;
    }

};