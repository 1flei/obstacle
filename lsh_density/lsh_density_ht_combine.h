#pragma once

#include "e2eigen.h"
#include "pert_hash.h"
#include <unordered_map>
#include <vector>
#include <fstream>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/unordered_map.hpp>

#include <functinoal>
#include "bottom_k_sketch.h"

#ifdef _MSC_VER
#  include <intrin.h>
#  define __builtin_popcount __popcnt
#endif

//using inclusion-exclusion principle
//using hash table which might be so efficient
class LSH_KDE
{
public:
    int dim, k, l, nhasher;
    E2Eigen lsh_hasher;
    std::vector<std::vector<std::unordered_map<E2Eigen::SigType, int32_t>> > hash_tabless;
    // will have l layers
    // the first layer has (m=nhasher) hashers
    // the second layer has m*(m-2)/2 hashers
    // the third layer hash m*(m-2)*(m-3)/6 hashers, and so on so forth


public:

    LSH_KDE(int dim, int k, int l, double r, int nhasher = 1)
        : dim(dim), k(k), l(l), nhasher(nhasher), lsh_hasher(dim, nhasher*k*l, r, k), hash_tabless(l), insert_fs(l)
    {
        assert(l < 10);
        //will not be efficient any more for larger l

        int64_t mchoicel = 1;
        for(int i=0;i<l;i++){
            mchoicel = mchoicel*(nhasher-i)/(i+1);
            hash_tabless[i].resize(mchoicel);
        }

    }

    //f:: layerid -> vecidx -> const vec<int>&  -> some IO
    template<class F>
    void for_m_choice_ls(int m, int l, std::vector<int>& curvecidx, 
        std::vector<int>& curvec, 
        const F& f)
    {
        if(curvec.size() < l){
            int beg = 0;
            if(!curvec.empty()) {
                f(curvecidx[curvec.size() -1], curvec);
                curvecidx[curvec.size() -1]++;

                beg = curvec.back()+1;
            }
            for(int i=beg; i< m-l; i++){
                curvec.push_back(i);
                for_m_choice_l(m, l, curvec);
                curvec.pop_back();
            }
        }
    }

    void insert(const std::vector<double> &x)
    {
        auto sig = lsh_hasher.getSigConcatenated(&x[0]);
        std::vector<int> curvecidx(l);
        std::fill(curvecidx.begin(), curvecidx.end(), 0);
        std::vector<int> curvec;

        auto f = [&](int vecidx, const std::vector<int>& vec){

        };

        //size of sig [nhasher*l]
        for(int i=0; i<nrepeat; i++) {
            auto& lsh_hasher = lsh_hashers[i];
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

                int hash_table_idx = i*l + j;
                hash_tables[hash_table_idx][hash_cache[mask]] += signs[mask];
            }
        }
    }

    double get_estimated_density(const std::vector<double> &x)
    {
        std::vector<double> cnts(nrepeat);
        for(int i=0; i<nrepeat; i++) {
            auto& lsh_hasher = lsh_hashers[i];
            auto sig = lsh_hasher.getSigConcatenated(&x[0]);
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

                int hash_table_idx = i*l + j;
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