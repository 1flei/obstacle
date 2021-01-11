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

#ifdef _MSC_VER
#  include <intrin.h>
#  define __builtin_popcount __popcnt
#endif

//using inclusion-exclusion principle
//using hash table which might be so efficient
class LSH_KDE_HT
{
public:
    int dim, k, l, nrepeat;
    std::vector<E2Eigen> lsh_hashers;
    std::vector<std::unordered_map<E2Eigen::SigType, int32_t>> hash_tables;
    // will have (2**l -1)*nrepeat hash tables

public:

    LSH_KDE_HT(int dim, int k, int l, double r, int nrepeat = 1)
        : dim(dim), k(k), l(l), nrepeat(nrepeat), hash_tables(((1<<l) -1)*nrepeat)
    {
        assert(l < 10);
        //will not be efficient any more for larger l
        lsh_hashers.reserve(nrepeat);
        for(int i=0;i<nrepeat;i++){
            lsh_hashers.emplace_back(dim, k*l, r, k, i*E2Eigen::DEFAULT_SEED);
        }
    }

    void insert(const std::vector<double> &x)
    {
        for(int i=0; i<nrepeat; i++) {
            auto& lsh_hasher = lsh_hashers[i];
            auto sig = lsh_hasher.get_sig_concatenated(&x[0]);
            std::vector<uint64_t> hash_cache(1<<l);
            hash_cache[0] = 0;
            std::vector<int> signs(1<<l);
            signs[0] = -1;

            //calculate hash_cache
            for(int j=0; j< (1<<l)-1; j++){
                uint32_t mask = j+1;
                uint32_t lsb = mask & -mask;
                uint32_t mask_prev = mask ^ lsb;
                int hashidx_to_combine = log2_1bit(lsb);
                hash_cache[mask] = hash_cache[mask_prev];
                hash_combine(hash_cache[mask], sig[hashidx_to_combine]);
                signs[mask] = signs[mask_prev] * -1;

                int hash_table_idx = i*((1<<l)-1) + j;
                hash_tables[hash_table_idx][hash_cache[mask]] += signs[mask];
            }
        }
    }

    double get_estimated_density(const std::vector<double> &x)
    {
        std::vector<double> cnts(nrepeat);
        for(int i=0; i<nrepeat; i++) {
            auto& lsh_hasher = lsh_hashers[i];
            auto sig = lsh_hasher.get_sig_concatenated(&x[0]);
            std::vector<uint64_t> hash_cache(1<<l);
            hash_cache[0] = 0;
            //calculate hash_cache
            for(int j=0; j< (1<<l)-1; j++){
                uint32_t mask = j+1;
                uint32_t lsb = mask & -mask;
                uint32_t mask_prev = mask ^ lsb;
                int hashidx_to_combine = log2_1bit(lsb);
                hash_cache[mask] = hash_cache[mask_prev];
                hash_combine(hash_cache[mask], sig[hashidx_to_combine]);

                int hash_table_idx = i*((1<<l)-1) + j;
                auto it = hash_tables[hash_table_idx].find(hash_cache[mask]);
                if(it != hash_tables[hash_table_idx].end()) {
                    cnts[i] += it->second;
                }
            }
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
        // std::vector<std::unordered_map<E2Eigen::SigType, int32_t>> hash_tables;
        for(int i=0;i<hash_tables.size();i++){
            ret += sizeof(hash_tables[i]);
            int64_t entry_size = sizeof(E2Eigen::SigType) + sizeof(int32_t) + sizeof(void*);
            int64_t bucket_size = sizeof(void*);
            ret += hash_tables[i].size()*entry_size + hash_tables[i].bucket_count() * bucket_size;
        }
        return ret;
    }


protected:

    int log2_1bit(uint32_t value)
    {
        static const int tab32[32] = {
            0, 9, 1, 10, 13, 21, 2, 29,
            11, 14, 16, 18, 22, 25, 3, 30,
            8, 12, 20, 28, 15, 17, 24, 7,
            19, 27, 23, 6, 26, 5, 4, 31};
        value = (value<<1) -1;
        return tab32[(uint32_t)(value * 0x07C4ACDD) >> 27];
    }
};