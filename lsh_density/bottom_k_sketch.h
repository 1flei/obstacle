#pragma once
#include "tabu_hash.h"
#include <boost/serialization/vector.hpp>
#include <boost/serialization/set.hpp>
#include <boost/container/flat_set.hpp>

//use boost::flat_set to maintain the bottom-k elems
namespace boost{
namespace serialization{
// flat_map
template<class Archive, class Type, class Key, class Compare, class Allocator >
inline void save(Archive & ar, const boost::container::flat_map<Key, Type, Compare, Allocator> &t, const unsigned int /* file_version */)
{
    boost::serialization::stl::save_collection<Archive,boost::container::flat_map<Key, Type, Compare, Allocator> >(ar, t);
}
template<class Archive, class Type, class Key, class Compare, class Allocator >
inline void load(Archive & ar, boost::container::flat_map<Key, Type, Compare, Allocator> &t, const unsigned int /* file_version */){
    load_map_collection(ar, t);
}
// split non-intrusive serialization function member into separate
// non intrusive save/load member functions
template<class Archive, class Type, class Key, class Compare, class Allocator >
inline void serialize(Archive & ar,boost::container::flat_map<Key, Type, Compare, Allocator> &t, const unsigned int file_version){
    boost::serialization::split_free(ar, t, file_version);
}
//flat_set
template<class Archive, class Key, class Compare, class Allocator >
inline void save(Archive & ar, const boost::container::flat_set<Key, Compare, Allocator> &t, const unsigned int /* file_version */){
    boost::serialization::stl::save_collection<Archive, boost::container::flat_set<Key, Compare, Allocator>>(ar, t);
}
template<class Archive, class Key, class Compare, class Allocator >
inline void load(Archive & ar, boost::container::flat_set<Key, Compare, Allocator> &t, const unsigned int /* file_version */){
    load_set_collection(ar, t);
}
// split non-intrusive serialization function member into separate
// non intrusive save/load member functions
template<class Archive, class Key, class Compare, class Allocator >
inline void serialize( Archive & ar, boost::container::flat_set<Key, Compare, Allocator> & t,const unsigned int file_version){
    boost::serialization::split_free(ar, t, file_version);
}
}}

template<class HASH_INT_TYPE>
class BottomKSketch;

template<class HASH_INT_TYPE>
class BottomKSketchSet
{
protected:
    // std::set<uint64_t> hashes;
    boost::container::flat_set<uint64_t> hashes;

    //if there is k elems, return maximum one
    //otherwise return 0
    uint64_t get_maximum_hash()
    {
        if (!hashes.empty()) 
            return *(hashes.rbegin());
        return 0;
    }

    void pop_maximum_hash()
    {
        auto it = std::prev(hashes.end());
        hashes.erase(it);
    }

    void insert(uint64_t hash)
    {
        hashes.insert(hash);
    }

    std::size_t size()
    {
        return hashes.size();
    }

    friend class BottomKSketch<HASH_INT_TYPE>;

public:
    template <class Archive>
    void serialize(Archive &ar, const unsigned int )
    {
        ar & hashes;
    }

    int64_t get_memory_usage()
    {
        return int64_t(sizeof(*this)) + hashes.capacity()*sizeof(uint64_t);
    }
};

template<class HASH_INT_TYPE>
class BottomKSketch
{
public:
    TabuHasher<HASH_INT_TYPE> hasher;
    int k;

    static const int64_t DEFAULT_SEED = 666;

    BottomKSketch(int k = 64, uint64_t seed = DEFAULT_SEED):
        hasher(seed), k(k)
    {
    }

    using Set = BottomKSketchSet<HASH_INT_TYPE>;

    void insert(HASH_INT_TYPE id, Set& hashes) 
    {
        uint64_t h = hasher.hash(id);
        if(hashes.size() < k || h < hashes.get_maximum_hash()){ 
            hashes.insert(h);
            if(hashes.size() > k) {
                hashes.pop_maximum_hash();
            }
        }
    }

    double get_cardinality(Set& hashes) 
    {
        if(hashes.size() < k){
            return hashes.size();
        } else {
            return (k-1)*1./hashes.get_maximum_hash()*TabuHasher<HASH_INT_TYPE>::max();
        }
    }
    
    Set union_set(const Set& a, const Set& b)
    {
        const Set& max_ab = a.size() > b.size() ? a : b;
        const Set& min_ab = a.size() > b.size() ? b : a;
        Set ret(max_ab);
        for(uint64_t k:min_ab) {
            ret.insert(k);
        }
        while(ret.size() > k){
            ret.pop_maximum_hash();
        }
        return ret;
    }
    
    //a = a cup b
    void merge(Set& a, const Set& b)
    {
        for(uint64_t k:b.hashes) {
            a.insert(k);
        }
        while(a.size() > k){
            a.pop_maximum_hash();
        }
    }

    int64_t get_memory_usage()
    {
        // TabuHasher<HASH_INT_TYPE> hasher;
        // int k;
        return int64_t(sizeof(int)) + hasher.get_memory_usage();
    }

public:
    template <class Archive>
    void serialize(Archive &ar, const unsigned int )
    {
        ar & hasher;
        ar & k;
    }
};

using BottomKSketch32 = BottomKSketch<int32_t>;