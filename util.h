#pragma once
#include <optional>
#include <vector>
#include "tl/function_ref.hpp"
#include <tuple>
#include <iomanip>

// const double time_interval = 600;  //ten minites
// const double time_interval_freq = 30;  //30 seconds
// const double time_interval = 30;  //30 seconds
// const int sliding_window_len = 6;
// const int sliding_window_step = 1;


const double eps = 1e-7;

inline bool eps_eq(double a, double b, double eps_=eps)
{
    return abs(a-b) < eps_;
}

inline bool eps_gte(double a, double b, double eps_=eps)
{
    return a-b > -eps_;
}

inline bool eps_gt(double a, double b, double eps_=eps)
{
    return a-b > eps_;
}


struct Pnt3
{
    double t, x, y;
    Pnt3() : t(0), x(0), y(0) {}
    Pnt3(double t, double x, double y): t(t), x(x), y(y) {}
    bool operator<(const Pnt3& p) const {
        return t<p.t;
    }
    template<class OS>
    friend OS& operator<<(OS& o, const Pnt3& p) {
        return o << "(" << p.t << "," << p.x << "," << p.y << ")";
    }
};

inline double eudist2(const Pnt3& p0, const Pnt3& p1)
{
    return (p0.x-p1.x)*(p0.x-p1.x) + (p0.y-p1.y)*(p0.y-p1.y);
}
inline double seg_eudist2(const Pnt3& p0, const Pnt3& p1, const Pnt3& p2, const Pnt3& p3)
{
    return (p1.x+p2.x-p0.x-p3.x)*(p1.x+p2.x-p0.x-p3.x) + (p1.y+p2.y-p0.y-p3.y)*(p1.y+p2.y-p0.y-p3.y);
}

struct Traj 
{
    std::vector<double> data;

    Pnt3& pnt_at(int i){
        return *((Pnt3*)&(data[i*3]));
    }

    Pnt3& operator[](int i){
        return pnt_at(i);
    }

    const Pnt3& pnt_at(int i) const{
        return *((Pnt3*)&(data[i*3]));
    }

    const Pnt3& operator[](int i) const{
        return pnt_at(i);
    }

    void push_back(const Pnt3& p){
        data.push_back(p.t);
        data.push_back(p.x);
        data.push_back(p.y);
    }
    void emplace_back(double t, double x, double y){
        data.push_back(t);
        data.push_back(x);
        data.push_back(y);
    }

    //might be unsafe
    void push_back_partial(double x){
        data.push_back(x);
    }

    std::vector<double>& get_flat_data(){
        return data;
    }

    int64_t get_num_pnts() const {
        return data.size()/3;
    }

    Pnt3* begin() const {
        return (Pnt3*)(&data[0]);
    }
    Pnt3* end() const {
        return begin()+get_num_pnts();
    }
};


struct TrajView
{
    Pnt3* p;
    int64_t nPnts;
    TrajView() : p(nullptr), nPnts(0) {}
    TrajView(const Pnt3* p_, int64_t nPnts_) : p((Pnt3*)p_), nPnts(nPnts_){}
    TrajView(const Traj& t, int offset, int64_t nPnts_) : p((Pnt3*)&t[offset]), nPnts(nPnts_){}

    Pnt3& pnt_at(int i){
        return *((Pnt3*)&(p[i]));
    }

    Pnt3& operator[](int i){
        return pnt_at(i);
    }

    const Pnt3& pnt_at(int i) const{
        return *((Pnt3*)&(p[i]));
    }

    const Pnt3& operator[](int i) const{
        return pnt_at(i);
    }

    int64_t get_num_pnts() const {
        return nPnts;
    }

    template<class OS>
    friend OS& operator<<(OS& o, const TrajView& idx) {
        if(idx.nPnts<=0) {
            return o << "[]";
        }

        o << "[" << idx.p[0];
        for(int i=1;i<idx.nPnts;i++){
            o << "," << idx.p[i];
        }
        return o << "]";
    }

    std::vector<double> to_flat_vec() const {
        std::vector<double> ret;
        for(int i=0;i<nPnts;i++){
            Pnt3& pt = p[i];
            ret.push_back(pt.t);
            ret.push_back(pt.x);
            ret.push_back(pt.y);
        }
        return ret;
    }
    std::vector<double> to_flat_vec_spatial() const {
        std::vector<double> ret;
        for(int i=0;i<nPnts;i++){
            Pnt3& pt = p[i];
            ret.push_back(pt.x);
            ret.push_back(pt.y);
        }
        return ret;
    }
    Pnt3* begin() const {
        return p;
    }
    Pnt3* end() const {
        return begin()+nPnts;
    }
};


inline std::optional<TrajView> try_get_trajview(const Traj& t, int offset, int64_t nPnts, double interval)
{
    if(offset <0 || offset+nPnts > t.get_num_pnts()) {
        return {};
    }
    double t0 = t[offset].t;
    double t1 = t[offset+nPnts-1].t;
    if( eps_eq(t1-t0, interval*(nPnts-1) ) ) {
        return TrajView(t, offset, nPnts);
    }
    return {};
}

inline std::optional<TrajView> try_get_trajview_clamp01(const Traj& t, int offset, int64_t nPnts, double interval)
{
    for(int cur_offset=offset; cur_offset>=0;--cur_offset){
        auto ret = try_get_trajview(t, cur_offset, nPnts, interval);
        if(ret){
            return ret;
        }
    }
    return {};
}

// Code from boost
// Reciprocal of the golden ratio helps spread entropy
//     and handles duplicates.
// See Mike Seymour in magic-numbers-in-boosthash-combine:
//     https://stackoverflow.com/questions/4948780

template <class T0, class T>
inline void hash_combine(T0& seed, T const& v)
{
    seed ^= std::hash<T>()(v) + 0x9e3779b9 + (seed<<6) + (seed>>2);
}

template <class T0, class T>
inline T0 hash_combine_bop(T0 const& seed, T const& v)
{
    return seed ^ (std::hash<T>()(v) + 0x9e3779b9 + (seed<<6) + (seed>>2) );
}

//---------------------------------------------------------
//tuple hasher
// function has to live in the std namespace 
// so that it is picked up by argument-dependent name lookup (ADL).
namespace std{
    namespace
    {
        // Recursive template code derived from Matthieu M.
        template <class Tuple, size_t Index = std::tuple_size<Tuple>::value - 1>
        struct HashValueImpl
        {
          static void apply(size_t& seed, Tuple const& tuple)
          {
            HashValueImpl<Tuple, Index-1>::apply(seed, tuple);
            hash_combine(seed, get<Index>(tuple));
          }
        };

        template <class Tuple>
        struct HashValueImpl<Tuple,0>
        {
          static void apply(size_t& seed, Tuple const& tuple)
          {
            hash_combine(seed, get<0>(tuple));
          }
        };
    }

    template <typename ... TT>
    struct hash<std::tuple<TT...>> 
    {
        size_t
        operator()(std::tuple<TT...> const& tt) const
        {                                              
            size_t seed = 0;                             
            HashValueImpl<std::tuple<TT...> >::apply(seed, tt);    
            return seed;                                 
        }                                              

    };
}

//format
#define FMT_HEADER_ONLY
#include "fmt/format.h"
#include "fmt/ranges.h"
#include "fmt/ostream.h"
#include "fmt/printf.h"

//use function_ref as callback
inline void sliding_window(const Traj& t,
    int len, 
    int step, 
    double interval, 
    const tl::function_ref<void(const TrajView&)>& f
)
{
    for(int i=0;i+len<=t.get_num_pnts();i+=step){
        double t0 = t[i].t;
        double t1 = t[i+len-1].t;
        if( eps_eq(t1-t0, interval*(len-1) ) ) {
            f(TrajView(&t[i], len));
        }
    }
}
inline void enum_sliding_window(const Traj& t,
    int len, 
    int step, 
    double interval, 
    const tl::function_ref<void(int, const TrajView&)>& f
)
{
    for(int i=0;i+len<=t.get_num_pnts();i+=step){
        double t0 = t[i].t;
        double t1 = t[i+len-1].t;
        // fmt::print("delta-t={}, actual_interval={}\n", t1-t0, interval*(len-1));
        if( eps_eq(t1-t0, interval*(len-1) ) ) {
            f(i, TrajView(&t[i], len));
        }
    }
}
inline void sliding_window_n(const Traj& t, 
    int cnt,
    int len, 
    int step, 
    double interval, 
    const tl::function_ref<void(const TrajView&)>& f
)
{
    for(int i=0;i+len<=t.get_num_pnts() && cnt>=0;i+=step){
        double t0 = t[i].t;
        double t1 = t[i+len-1].t;
        if( eps_eq(t1-t0, interval*(len-1) ) ) {
            f(TrajView(&t[i], len));
            --cnt;
        }
    }
}

// -----------------------------------------------------------------------------
template<class ScalarType>
ScalarType sqr(ScalarType x)
{
	return x*x;
}

//FProd :: ParamType -> ParamType -> ScalarType
//FSum :: ScalarType -> ScalarType -> ScalarType
template<class ScalarType, class ParamType, class FProd, class FSum> 
inline ScalarType fast_reduce(int dim, const ParamType* x, const ParamType* y, 
        FProd&& fp, FSum&& fs)
{
    unsigned d = dim & ~unsigned(7);
    const ParamType *aa = x, *end_a = aa + d;
    const ParamType *bb = y, *end_b = bb + d;
#ifdef __GNUC__
    __builtin_prefetch(aa, 0, 3);
    __builtin_prefetch(bb, 0, 0);
#endif
    ScalarType r = 0.0;
    ScalarType r0, r1, r2, r3, r4, r5, r6, r7;

    const ParamType *a = end_a, *b = end_b;

    r0 = r1 = r2 = r3 = r4 = r5 = r6 = r7 = 0.0;

    switch (dim & 7) {
        case 7: r6 = fp(a[6], b[6]);
  		// fall through
        case 6: r5 = fp(a[5], b[5]);
  		// fall through
        case 5: r4 = fp(a[4], b[4]);
  		// fall through
        case 4: r3 = fp(a[3], b[3]);
  		// fall through
        case 3: r2 = fp(a[2], b[2]);
  		// fall through
        case 2: r1 = fp(a[1], b[1]);
  		// fall through
        case 1: r0 = fp(a[0], b[0]);
    }

    a = aa; b = bb;
	const auto fsum8 = [&](){
		auto r01 = fs(r0, r1);
		auto r23 = fs(r2, r3);
		auto r45 = fs(r4, r5);
		auto r67 = fs(r6, r7);
		auto r0123 = fs(r01, r23);
		auto r4567 = fs(r45, r67);
		return fs(r0123, r4567);
	};

    for (; a < end_a; a += 8, b += 8) {
#ifdef __GNUC__
        __builtin_prefetch(a + 32, 0, 3);
        __builtin_prefetch(b + 32, 0, 0);
#endif
		r = fs(r, fsum8() );
		r0 = fp(a[0], b[0]);
		r1 = fp(a[1], b[1]);
		r2 = fp(a[2], b[2]);
		r3 = fp(a[3], b[3]);
		r4 = fp(a[4], b[4]);
		r5 = fp(a[5], b[5]);
		r6 = fp(a[6], b[6]);
		r7 = fp(a[7], b[7]);
    }

    r = fs(r, fsum8() );
    return r;
}


template<class ScalarType> 
inline ScalarType calc_l2_sqr(int dim, const ScalarType* x, const ScalarType* y)
{
	const auto fProd = [](ScalarType a, ScalarType b) -> ScalarType{
		return sqr(a-b);
	};
	const auto fSum = [](ScalarType a, ScalarType b) -> ScalarType{
		return a+b;
	};
	return fast_reduce<ScalarType, ScalarType>(dim, x, y, fProd, fSum);
}


template<class ScalarType> 
inline ScalarType calc_l1_dist(int dim, const ScalarType* x, const ScalarType* y)
{
	const auto fProd = [](ScalarType a, ScalarType b){
		return abs(a-b);
	};
	const auto fSum = [](ScalarType a, ScalarType b){
		return a+b;
	};
	return fast_reduce<ScalarType, ScalarType>(dim, x, y, fProd, fSum);
}

template<class ScalarType> 
inline ScalarType calc_inner_product(int dim, const ScalarType* x, const ScalarType* y)
{
	const auto fProd = [](ScalarType a, ScalarType b){
		return a*b;
	};
	const auto fSum = [](ScalarType a, ScalarType b){
		return a+b;
	};
	return fast_reduce<ScalarType, ScalarType>(dim, x, y, fProd, fSum);
}

// -----------------------------------------------------------------------------
template<class ScalarType>
ScalarType calc_l2_dist(					// calc L2 distance
	int   dim,							// dimension
	const ScalarType *p1,					// 1st point
	const ScalarType *p2)					// 2nd point
{
	return sqrt(calc_l2_sqr(dim, p1, p2));
}


struct BBox
{
    double minlat, maxlat;
    double minlng, maxlng;
    double mint, maxt;


    BBox() : minlat(-180), maxlat(180), minlng(-180), maxlng(180), mint(0), maxt(2e9){}
    BBox(double minlng, double maxlng, double minlat, double maxlat, double mint, double maxt) : 
        minlat(minlat), maxlat(maxlat), minlng(minlng), maxlng(maxlng), mint(mint), maxt(maxt){}
};



template<class IdxType=int> 
class HashDisjointSet
{
public:
	HashDisjointSet(){}
	~HashDisjointSet(){}

	IdxType get(IdxType u){
		if(parent.find(u)==parent.end()){
			parent.emplace(u, u);
			return u;
		} else if(parent[u]==u){
			return u;
		} else{
			IdxType ret = get(parent[u]);
			parent[u] = ret;
			return ret;
		}
	}

	//make sure pu is parent of pv
	void merge(IdxType u, IdxType v){
		IdxType pu = get(u);
		IdxType pv = get(v);
		if(pu!=pv){
			parent[pv] = pu;
		}
	}
	std::unordered_map<IdxType, IdxType> parent;
};