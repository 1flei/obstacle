#pragma once
#include <string> 
#include <memory>
#include "hnswlib/hnswlib.h"
#include "hnswlib/hnswalg.h"
#include "subtidx.h"
#include "my_pri_queue.h"
#include <thread>

inline double PntEuDist(const void *pVect1v, const void *pVect2v, const void *qty_ptr) 
{
    Pnt3 *pVect1 = (Pnt3 *) pVect1v;
    Pnt3 *pVect2 = (Pnt3 *) pVect2v;
    size_t dim = *((size_t *) qty_ptr);
    auto fsum = [](double a, double b) -> double {
        return a+b;
    };

    return fast_reduce<double, Pnt3>(dim, pVect1, pVect2, eudist2, fsum);
}

class Pnt3EuSpace : public hnswlib::SpaceInterface<double> 
{
    size_t data_size;
    size_t dim;
public:
    Pnt3EuSpace(size_t dim) :dim(dim) {
        data_size = sizeof(Pnt3) * dim;
    }

    size_t get_data_size() {
        return data_size;
    }

    hnswlib::DISTFUNC<double> get_dist_func() {
        return PntEuDist;
    }

    void *get_dist_func_param() {
        return &dim;
    }

    ~Pnt3EuSpace() {}
};

const int band_width = 2;
inline double PntDTW(const void *pVect1v, const void *pVect2v, const void *qty_ptr) 
{
    Pnt3 *pVect1 = (Pnt3 *) pVect1v;
    Pnt3 *pVect2 = (Pnt3 *) pVect2v;
    size_t dim = *((size_t *) qty_ptr);

    static std::vector<double> dtws((dim+1) * (dim+1));
    dtws.resize((dim+1) * (dim+1));
    auto dtwat= [&](int i, int j)->double&{
        return dtws[i*dim + j];
    };
    const double inf = 1e9;
    for(int i=0;i<dim+1;i++){
        dtwat(std::max(0, i-band_width-1), i) = inf;
        dtwat(i, std::max(0, i-band_width-1)) = inf;
    }
    dtwat(0, 0) = 0.;

    for(int i=1;i<dim+1;i++){
        for(int j=std::max(1, i-band_width); j < std::min<int>(dim, i+band_width)+1; j++) {
            double cost = eudist2(pVect1[i-1], pVect2[j-1]);
            dtwat(i, j) = cost + std::min(dtwat(i-1, j), std::min(dtwat(i, j-1), dtwat(i-1, j-1)));
        }
    }
    return dtwat(dim, dim);
}
inline double DiscreteFrechet(const void *pVect1v, const void *pVect2v, const void *qty_ptr) 
{
    Pnt3 *pVect1 = (Pnt3 *) pVect1v;
    Pnt3 *pVect2 = (Pnt3 *) pVect2v;
    size_t dim = *((size_t *) qty_ptr);

    static std::vector<double> frechets((dim+1) * (dim+1));
    frechets.resize((dim+1) * (dim+1));
    auto distat= [&](int i, int j)->double&{
        return frechets[i*dim + j];
    };
    const double inf = 1e9;
    for(int i=0;i<dim+1;i++){
        distat(std::max(0, i-band_width-1), i) = inf;
        distat(i, std::max(0, i-band_width-1)) = inf;
    }
    distat(0, 0) = 0.;

    for(int i=1;i<dim+1;i++){
        for(int j=std::max(1, i-band_width); j < std::min<int>(dim, i+band_width)+1; j++) {
            double cost = eudist2(pVect1[i-1], pVect2[j-1]);
            distat(i, j) = std::max(cost, std::min(distat(i-1, j), std::min(distat(i, j-1), distat(i-1, j-1))) );
        }
    }
    return distat(dim, dim);
}
inline double DiscreteFrechetNormalized(const void *pVect1v, const void *pVect2v, const void *qty_ptr) 
{
    Pnt3 *pVect1 = (Pnt3 *) pVect1v;
    Pnt3 *pVect2 = (Pnt3 *) pVect2v;
    size_t dim = *((size_t *) qty_ptr);
    double dtw = DiscreteFrechet(pVect1v, pVect2v, qty_ptr);

    const double eps = 1e-9;
    double lenp1 = eps;
    double lenp2 = eps;
    for(int i=0;i<dim-1;i++){
        lenp1 += sqrt(eudist2(pVect1[i], pVect1[i+1]));
        lenp2 += sqrt(eudist2(pVect2[i], pVect2[i+1]));
    }
    return dtw/lenp1/lenp2;
}

inline double PntDTWNormalized(const void *pVect1v, const void *pVect2v, const void *qty_ptr) 
{
    Pnt3 *pVect1 = (Pnt3 *) pVect1v;
    Pnt3 *pVect2 = (Pnt3 *) pVect2v;
    size_t dim = *((size_t *) qty_ptr);
    double dtw = PntDTW(pVect1v, pVect2v, qty_ptr);

    const double eps = 1e-9;
    double lenp1 = eps;
    double lenp2 = eps;
    for(int i=0;i<dim-1;i++){
        lenp1 += sqrt(eudist2(pVect1[i], pVect1[i+1]));
        lenp2 += sqrt(eudist2(pVect2[i], pVect2[i+1]));
    }
    return dtw/lenp1/lenp2;
}


class Pnt3DTWSpace : public hnswlib::SpaceInterface<double> 
{
    size_t data_size;
    size_t dim;
public:
    Pnt3DTWSpace(size_t dim) :dim(dim) {
        data_size = sizeof(Pnt3) * dim;
    }

    size_t get_data_size() {
        return data_size;
    }

    hnswlib::DISTFUNC<double> get_dist_func() {
        return PntDTW;
    }

    void *get_dist_func_param() {
        return &dim;
    }

    ~Pnt3DTWSpace() {}
};

class Pnt3DTWNormalizedSpace : public hnswlib::SpaceInterface<double> 
{
    size_t data_size;
    size_t dim;
public:
    Pnt3DTWNormalizedSpace(size_t dim) :dim(dim) {
        data_size = sizeof(Pnt3) * dim;
    }

    size_t get_data_size() {
        return data_size;
    }

    hnswlib::DISTFUNC<double> get_dist_func() {
        return PntDTWNormalized;
    }

    void *get_dist_func_param() {
        return &dim;
    }

    ~Pnt3DTWNormalizedSpace() {}
};

class DiscreteFrechetSpace : public hnswlib::SpaceInterface<double> 
{
    size_t data_size;
    size_t dim;
public:
    DiscreteFrechetSpace(size_t dim) :dim(dim) {
        data_size = sizeof(Pnt3) * dim;
    }

    size_t get_data_size() {
        return data_size;
    }

    hnswlib::DISTFUNC<double> get_dist_func() {
        return DiscreteFrechet;
    }

    void *get_dist_func_param() {
        return &dim;
    }

    ~DiscreteFrechetSpace() {}
};

class DiscreteFrechetNormalizedSpace : public hnswlib::SpaceInterface<double> 
{
    size_t data_size;
    size_t dim;
public:
    DiscreteFrechetNormalizedSpace(size_t dim) :dim(dim) {
        data_size = sizeof(Pnt3) * dim;
    }

    size_t get_data_size() {
        return data_size;
    }

    hnswlib::DISTFUNC<double> get_dist_func() {
        return DiscreteFrechet;
    }

    void *get_dist_func_param() {
        return &dim;
    }

    ~DiscreteFrechetNormalizedSpace() {}
};

void build_index(std::unordered_map<int, Traj>& traj_dict, 
    int len, 
    int step, 
    double interval, 
    hnswlib::HierarchicalNSW<double, SubtIdx>& hnsw, 
    int n_threads=16);

SubtPriQueue searchKNN(hnswlib::HierarchicalNSW<double, SubtIdx>& hnsw, void* query, int k);