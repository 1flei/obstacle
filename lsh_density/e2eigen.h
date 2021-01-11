#pragma once

#include <vector>
#include <cassert>
#include <random>
#include "../Eigen/Eigen"
#include "../util.h"
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/unordered_map.hpp>


namespace boost {
   namespace serialization {
      template <class Archive> 
      void serialize( Archive & ar, Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& g, const unsigned int ){
          ar & boost::serialization::make_array(g.derived().data(), g.size());
      }
      template <class Archive> 
      void serialize( Archive & ar, Eigen::Matrix<double, Eigen::Dynamic, 1>& g, const unsigned int ){
          ar & boost::serialization::make_array(g.derived().data(), g.size());
      }
   }
}

//Random Projection using eigen
class E2Eigen
{
public:
    using SigType = int32_t;
    static const int64_t DEFAULT_SEED = 666;
    E2Eigen()
    {
        //for serialization only, 
        //TODO: will be modified to save_construct and load_construct in the future
    }
    E2Eigen(int d, int K, double r, int nComb=1, int64_t seed = DEFAULT_SEED)      //dim of data object, #hasher, radius 
        :dim(d), K(K), r(r), nComb(nComb)
    {
        assert(d > 0 && K > 0);

        std::normal_distribution<double> normal(0.);
        std::uniform_real_distribution<double> uniform(0., r);
        // std::random_device rd;
        // std::default_random_engine rng(rd());
        std::default_random_engine rng(seed);

        p.resize(K, d);
        b.resize(K);
        for (int i = 0; i < K; i++) {
            for(int j=0;j<d;j++){
                p(i, j) = normal(rng);
            }
        }
        for (int i = 0; i < K; i++) {
            b(i) = uniform(rng) + 0.5*r;
        }
        sigdim = K;
    }
    ~E2Eigen() {}
    std::vector<SigType> get_sig(const double *data) const 
    {
        std::vector<SigType> ret(sigdim);
        get_sig(data, &ret[0]);
        return ret;
    }
    void get_sig(const double *data, SigType* ret) const 
    {
        //assume memory is allocated for ret
        Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, 1> > data_v(data, dim);
        Eigen::Map<Eigen::Matrix<SigType, Eigen::Dynamic, 1> > ret_v(ret, sigdim);

        auto ret_r = (p*data_v+b)/r;
        ret_v = ret_r.cast<int32_t>();
    }

    std::vector<SigType> get_sig_concatenated(const double *data) const 
    {
        std::vector<SigType> sigraw = get_sig(data);

        std::vector<SigType> sig(sigdim/ nComb);
        for(int i=0;i<sigdim / nComb; i++){
            uint64_t sigi = 0;
            for(int j=0; j<nComb; j++){
                int  sigidx = i*nComb + j;
                hash_combine(sigi, sigraw[sigidx]);
            }
            sig[i] = sigi;
        }
        return sig;
    }

    int64_t get_memory_usage()
    {
        return int64_t(sizeof(*this)) + (K*dim)*sizeof(double) + int64_t(K)*sizeof(double);
    }

    template <class Archive>
    void serialize(Archive &ar, const unsigned int )
    {
        ar & dim;
        ar & K;
        ar & r;
        ar & sigdim;
        ar & nComb;
        ar & p;
        ar & b;
    }

    int dim, K;
    double r;
    int sigdim;
    int nComb;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> p;
    Eigen::Matrix<double, Eigen::Dynamic, 1> b;
    
    // std::vector<double> p;
    // std::vector<double> b;
};


