#pragma once 


#include "../lsh_density/lsh_density.h"

class LSH_KDE_PATCH
{
public:
    LSH_KDE* lsh_kde;

    LSH_KDE_PATCH(LSH_KDE& lsh_kde)
        :lsh_kde(&lsh_kde), incremental_hash_tables(lsh_kde.hash_tables.size())
    {} 

    // std::vector<std::vector<int32_t>> hash_tables;
    std::vector<std::unordered_map<int32_t, int32_t> > incremental_hash_tables;

    void insert(const std::vector<double> &x, double inc=1.)
    {
        // fmt::print("x={}\n", x);
        for(int i=0; i<lsh_kde->nrepeat; i++) {
            auto& lsh_hasher = lsh_kde->lsh_hashers[i];
            auto sig = lsh_hasher.get_sig_concatenated(&x[0]);

            for(int k=0;k<lsh_kde->count_min_sketch_nrepeat;k++){
                std::vector<uint64_t> hash_cache(1<<lsh_kde->l);
                hash_cache[0] = 0;

                //calculate hash_cache
                for(int j=0; j< lsh_kde->expl; j++){
                    uint32_t mask = j+1;
                    uint32_t lsb = mask & -mask;
                    uint32_t mask_prev = mask ^ lsb;
                    int hashidx_to_combine = lsh_kde->log2_1bit(lsb);
                    // hash_cache[mask] = hash_cache[mask_prev];
                    // hash_combine(hash_cache[mask], sig[hashidx_to_combine]);
                    hash_cache[mask] = hash_cache[mask_prev] ^ lsh_kde->pert_hashers[k*lsh_kde->l+hashidx_to_combine].hash(sig[hashidx_to_combine]);

                    int hash_table_idx = i*lsh_kde->expl + j;
                    incremental_hash_tables[hash_table_idx][hash_cache[mask]%lsh_kde->ht_size] += inc;
                }
            }
        }
    }

    void insert_if_not_exist(const std::vector<double> &x, double inc= 1.)
    {
        // fmt::print("x={}\n", x);
        for(int i=0; i<lsh_kde->nrepeat; i++) {
            auto& lsh_hasher = lsh_kde->lsh_hashers[i];
            auto sig = lsh_hasher.get_sig_concatenated(&x[0]);

            // fmt::print("sig={}\n", sig);

            for(int k=0;k<lsh_kde->count_min_sketch_nrepeat;k++){
                std::vector<uint64_t> hash_cache(1<<lsh_kde->l);
                hash_cache[0] = 0;

                //calculate hash_cache
                for(int j=0; j< lsh_kde->expl; j++){
                    uint32_t mask = j+1;
                    uint32_t lsb = mask & -mask;
                    uint32_t mask_prev = mask ^ lsb;
                    int hashidx_to_combine = lsh_kde->log2_1bit(lsb);
                    // hash_cache[mask] = hash_cache[mask_prev];
                    // hash_combine(hash_cache[mask], sig[hashidx_to_combine]);
                    hash_cache[mask] = hash_cache[mask_prev] ^ lsh_kde->pert_hashers[k*lsh_kde->l+hashidx_to_combine].hash(sig[hashidx_to_combine]);

                    int hash_table_idx = i*lsh_kde->expl + j;
                    int key = hash_cache[mask]%lsh_kde->ht_size;
                    if(incremental_hash_tables[hash_table_idx].find(key) == incremental_hash_tables[hash_table_idx].end()) {
                        incremental_hash_tables[hash_table_idx][key] += inc;
                    }
                }
            }
        }
    }

    void merge_to_lsh_kde() {
        for(int i=0;i<incremental_hash_tables.size();i++){
            auto& lsh_ht = lsh_kde->hash_tables[i];
            auto& inc_ht = incremental_hash_tables[i];
            for(auto&& p:inc_ht){
                lsh_ht[p.first] += p.second;
            }
        }
    }
};



class LSH_KDE_SUCC_PATCH
{
public:
    LSH_KDE* lsh_kde;

    LSH_KDE_SUCC_PATCH(LSH_KDE& lsh_kde)
        :lsh_kde(&lsh_kde), incremental_hash_tables(lsh_kde.hash_tables.size())
    {} 

    // std::vector<std::vector<int32_t>> hash_tables;
    std::vector<std::unordered_map<int32_t, int32_t> > incremental_hash_tables;

    void insert(const std::vector<double> &x, const std::vector<double> &y, double inc=1.)
    {
        // fmt::print("x={}\n", x);
        for(int i=0; i<lsh_kde->nrepeat; i++) {
            //use different hash for x and y
            auto& lsh_hasherx = lsh_kde->lsh_hashers[i];
            auto& lsh_hashery = lsh_kde->lsh_hashers[(i+1)%lsh_kde->nrepeat];
            auto sigx = lsh_hasherx.get_sig_concatenated(&x[0]);
            auto sigy = lsh_hashery.get_sig_concatenated(&y[0]);

            for(int k=0;k<lsh_kde->count_min_sketch_nrepeat;k++){
                std::vector<uint64_t> hash_cache(1<<lsh_kde->l);
                hash_cache[0] = 0;

                //calculate hash_cache
                for(int j=0; j< lsh_kde->expl; j++){
                    uint32_t mask = j+1;
                    uint32_t lsb = mask & -mask;
                    uint32_t mask_prev = mask ^ lsb;
                    int hashidx_to_combine = lsh_kde->log2_1bit(lsb);
                    // hash_cache[mask] = hash_cache[mask_prev];
                    // hash_combine(hash_cache[mask], sigx[hashidx_to_combine], sigy[hashidx_to_combine]);
                    auto& cur_pert_hasher = lsh_kde->pert_hashers[k*lsh_kde->l+hashidx_to_combine];
                    hash_cache[mask] = hash_cache[mask_prev] ^ cur_pert_hasher.hash(sigx[hashidx_to_combine]) ^
                            cur_pert_hasher.hash(sigy[hashidx_to_combine]);

                    int hash_table_idx = i*lsh_kde->expl + j;
                    incremental_hash_tables[hash_table_idx][hash_cache[mask]%lsh_kde->ht_size] += inc;
                }
            }
        }
    }

    void insert_if_not_exist(const std::vector<double> &x, const std::vector<double> &y, double inc= 1.)
    {
        // fmt::print("x={}\n", x);
        for(int i=0; i<lsh_kde->nrepeat; i++) {
            //use different hash for x and y
            auto& lsh_hasherx = lsh_kde->lsh_hashers[i];
            auto& lsh_hashery = lsh_kde->lsh_hashers[(i+1)%lsh_kde->nrepeat];
            auto sigx = lsh_hasherx.get_sig_concatenated(&x[0]);
            auto sigy = lsh_hashery.get_sig_concatenated(&y[0]);

            for(int k=0;k<lsh_kde->count_min_sketch_nrepeat;k++){
                std::vector<uint64_t> hash_cache(1<<lsh_kde->l);
                hash_cache[0] = 0;

                //calculate hash_cache
                for(int j=0; j< lsh_kde->expl; j++){
                    uint32_t mask = j+1;
                    uint32_t lsb = mask & -mask;
                    uint32_t mask_prev = mask ^ lsb;
                    int hashidx_to_combine = lsh_kde->log2_1bit(lsb);
                    // hash_cache[mask] = hash_cache[mask_prev];
                    // hash_combine(hash_cache[mask], sigx[hashidx_to_combine], sigy[hashidx_to_combine]);
                    auto& cur_pert_hasher = lsh_kde->pert_hashers[k*lsh_kde->l+hashidx_to_combine];
                    hash_cache[mask] = hash_cache[mask_prev] ^ cur_pert_hasher.hash(sigx[hashidx_to_combine]) ^
                            cur_pert_hasher.hash(sigy[hashidx_to_combine]);

                    int hash_table_idx = i*lsh_kde->expl + j;
                    int key = hash_cache[mask]%lsh_kde->ht_size;
                    if(incremental_hash_tables[hash_table_idx].find(key) == incremental_hash_tables[hash_table_idx].end()) {
                        incremental_hash_tables[hash_table_idx][key] += inc;
                    }
                }
            }
        }
    }

    void merge_to_lsh_kde() {
        for(int i=0;i<incremental_hash_tables.size();i++){
            auto& lsh_ht = lsh_kde->hash_tables[i];
            auto& inc_ht = incremental_hash_tables[i];
            for(auto&& p:inc_ht){
                lsh_ht[p.first] += p.second;
            }
        }
    }
};


inline void insert_succ_density(LSH_KDE& lsh_kde, const std::vector<double> &x, const std::vector<double> &y, double inc=1.)
{
    // fmt::print("x={}\n", x);
    for(int i=0; i<lsh_kde.nrepeat; i++) {
        //use different hash for x and y
        auto& lsh_hasherx = lsh_kde.lsh_hashers[i];
        auto& lsh_hashery = lsh_kde.lsh_hashers[(i+1)%lsh_kde.nrepeat];
        auto sigx = lsh_hasherx.get_sig_concatenated(&x[0]);
        auto sigy = lsh_hashery.get_sig_concatenated(&y[0]);

        for(int k=0;k<lsh_kde.count_min_sketch_nrepeat;k++){
            std::vector<uint64_t> hash_cache(1<<lsh_kde.l);
            hash_cache[0] = 0;

            //calculate hash_cache
            for(int j=0; j< lsh_kde.expl; j++){
                uint32_t mask = j+1;
                uint32_t lsb = mask & -mask;
                uint32_t mask_prev = mask ^ lsb;
                int hashidx_to_combine = lsh_kde.log2_1bit(lsb);
                // hash_cache[mask] = hash_cache[mask_prev];
                // hash_combine(hash_cache[mask], sigx[hashidx_to_combine], sigy[hashidx_to_combine]);
                auto& cur_pert_hasher = lsh_kde.pert_hashers[k*lsh_kde.l+hashidx_to_combine];
                hash_cache[mask] = hash_cache[mask_prev] ^ cur_pert_hasher.hash(sigx[hashidx_to_combine]) ^
                        cur_pert_hasher.hash(sigy[hashidx_to_combine]);

                int hash_table_idx = i*lsh_kde.expl + j;
                lsh_kde.hash_tables[hash_table_idx][hash_cache[mask]%lsh_kde.ht_size] += inc;
            }
        }
    }
}

inline double get_estimated_succ_density(const LSH_KDE& lsh_kde, const std::vector<double>& x, const std::vector<double>& y)
{
    int count_min_sketch_nrepeat = lsh_kde.count_min_sketch_nrepeat;
    int nrepeat = lsh_kde.nrepeat;
    int expl = lsh_kde.expl;
    int l = lsh_kde.l;
    int ht_size = lsh_kde.ht_size;

    std::vector<double> cnts(nrepeat*(expl)*count_min_sketch_nrepeat);

    for(int i=0; i<nrepeat; i++) {
        auto& lsh_hasherx = lsh_kde.lsh_hashers[i];
        auto& lsh_hashery = lsh_kde.lsh_hashers[(i+1)%lsh_kde.nrepeat];
        auto sigx = lsh_hasherx.get_sig_concatenated(&x[0]);
        auto sigy = lsh_hashery.get_sig_concatenated(&y[0]);
        // fmt::print("for hasher-{}, sig={}\n", i, sig);
        for(int k=0;k<count_min_sketch_nrepeat;k++){
            std::vector<uint64_t> hash_cache(1<<l);
            hash_cache[0] = 0;
            //calculate hash_cache
            for(int j=0; j< expl; j++){
                uint32_t mask = j+1;
                uint32_t lsb = mask & -mask;
                uint32_t mask_prev = mask ^ lsb;
                int hashidx_to_combine = lsh_kde.log2_1bit(lsb);                    
                // hash_cache[mask] = hash_cache[mask_prev];
                // hash_combine(hash_cache[mask], sig[hashidx_to_combine]);
                auto& cur_pert_hasher = lsh_kde.pert_hashers[k*lsh_kde.l+hashidx_to_combine];
                hash_cache[mask] = hash_cache[mask_prev] ^ cur_pert_hasher.hash(sigx[hashidx_to_combine]) ^
                        cur_pert_hasher.hash(sigy[hashidx_to_combine]);

                int hash_table_idx = i*expl + j;
                // cnts[i] += hash_tables[hash_table_idx][hash_cache[mask]%ht_size];

                cnts[(i*(expl)+j)*count_min_sketch_nrepeat+k] = lsh_kde.hash_tables[hash_table_idx][hash_cache[mask]%ht_size];
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
            current_density_estimator += lsh_kde.hash_signs[j+1] * mincnt;
        }
        sum += current_density_estimator;
    }
    return sum/nrepeat;
}