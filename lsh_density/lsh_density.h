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
class LSH_KDE
{
public:
    int dim, k, l, nrepeat, ht_size, count_min_sketch_nrepeat, expl;
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

    std::vector<int> hash_signs;

public:

    LSH_KDE(int dim, int k, int l, double r, int nrepeat = 1, int ht_size = 4194301, int count_min_sketch_nrepeat = 1)
        : dim(dim), k(k), l(l), nrepeat(nrepeat), ht_size(ht_size), count_min_sketch_nrepeat(count_min_sketch_nrepeat), expl((1<<l)-1), hash_tables(((1<<l) -1)*nrepeat)
    {
        assert(l < 10);
        //will not be efficient any more for larger l
        lsh_hashers.reserve(nrepeat);
        for(int i=0;i<nrepeat;i++){
            lsh_hashers.emplace_back(dim, k*l, r, k, i*E2Eigen::DEFAULT_SEED);
        }

        for(int i=0;i<hash_tables.size();i++){
            hash_tables[i].resize(ht_size);
        }

        pert_hashers.reserve(count_min_sketch_nrepeat*l);
        for(int i=0;i<count_min_sketch_nrepeat*l;i++){
            pert_hashers.emplace_back(i*TabuHasher32::DEFAULT_SEED);
        }

        //compute signs
        hash_signs.resize(1<<l);
        hash_signs[0] = -1;
        for(int j=0;j<expl; j++){
            uint32_t mask = j+1;
            uint32_t lsb = mask & -mask;
            uint32_t mask_prev = mask ^ lsb;
            hash_signs[mask] = hash_signs[mask_prev] * -1;
        }
    }

    void insert(const std::vector<double> &x)
    {
        // fmt::print("x={}\n", x);
        for(int i=0; i<nrepeat; i++) {
            auto& lsh_hasher = lsh_hashers[i];
            auto sig = lsh_hasher.get_sig_concatenated(&x[0]);

            for(int k=0;k<count_min_sketch_nrepeat;k++){
                std::vector<uint64_t> hash_cache(1<<l);
                hash_cache[0] = 0;

                //calculate hash_cache
                for(int j=0; j< expl; j++){
                    uint32_t mask = j+1;
                    uint32_t lsb = mask & -mask;
                    uint32_t mask_prev = mask ^ lsb;
                    int hashidx_to_combine = log2_1bit(lsb);
                    // hash_cache[mask] = hash_cache[mask_prev];
                    // hash_combine(hash_cache[mask], sig[hashidx_to_combine]);
                    hash_cache[mask] = hash_cache[mask_prev] ^ pert_hashers[k*l+hashidx_to_combine].hash(sig[hashidx_to_combine]);

                    int hash_table_idx = i*expl + j;
                    hash_tables[hash_table_idx][hash_cache[mask]%ht_size] += 1;
                }
            }
        }
    }

    double get_estimated_density(const std::vector<double> &x)
    {
        // fmt::print("x={}\n", x);
        std::vector<double> cnts(nrepeat*(expl)*count_min_sketch_nrepeat);

        for(int i=0; i<nrepeat; i++) {
            auto& lsh_hasher = lsh_hashers[i];
            auto sig = lsh_hasher.get_sig_concatenated(&x[0]);
            // fmt::print("for hasher-{}, sig={}\n", i, sig);
            for(int k=0;k<count_min_sketch_nrepeat;k++){
                std::vector<uint64_t> hash_cache(1<<l);
                hash_cache[0] = 0;
                //calculate hash_cache
                for(int j=0; j< expl; j++){
                    uint32_t mask = j+1;
                    uint32_t lsb = mask & -mask;
                    uint32_t mask_prev = mask ^ lsb;
                    int hashidx_to_combine = log2_1bit(lsb);                    
                    // hash_cache[mask] = hash_cache[mask_prev];
                    // hash_combine(hash_cache[mask], sig[hashidx_to_combine]);
                    hash_cache[mask] = hash_cache[mask_prev] ^ pert_hashers[k*l+hashidx_to_combine].hash(sig[hashidx_to_combine]);

                    int hash_table_idx = i*expl + j;
                    // cnts[i] += hash_tables[hash_table_idx][hash_cache[mask]%ht_size];

                    cnts[(i*(expl)+j)*count_min_sketch_nrepeat+k] = hash_tables[hash_table_idx][hash_cache[mask]%ht_size];
                }
            }
        }

        //mean over nrepeat
        double sum = 0;
        for(int i=0;i<nrepeat;i++){
            //sum the cnt up using the exclusion and inclusion principle
            double current_density_estimator = 0;
            for(int j=0;j<expl;j++) {
                //count-min sketch like estimator
                int mincnt = 1e8;
                for(int k=0;k<count_min_sketch_nrepeat;k++){
                    int idx = (i*(expl)+j)*count_min_sketch_nrepeat+k;
                    if(cnts[idx] < mincnt){
                        mincnt = cnts[idx];
                    }
                    // printf("i=%d, j=%d, k=%d, cnt=%f\n", i, j, k, cnts[idx]);
                }
                current_density_estimator += hash_signs[j+1] * mincnt;
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
        ar & k;
        ar & l;
        ar & nrepeat;
        ar & ht_size;
        ar & count_min_sketch_nrepeat;
        ar & expl;
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
        // std::vector<int> hash_signs;
        ret += hash_signs.capacity() *sizeof(int);
        return ret;
    }

    int log2_1bit(uint32_t value) const 
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

