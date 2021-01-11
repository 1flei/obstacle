#pragma once 

#include <random>
#include <array>
#include <cinttypes>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/array.hpp>


//combine hash signature using tabulation hash
template<class HASH_INT_TYPE>
struct TabuHasher
{
    static const int64_t DEFAULT_SEED = 666;
	static const int SIZEOF_HASH_INT_TYPE = sizeof(HASH_INT_TYPE);

	std::array<std::array<uint64_t, 256>, SIZEOF_HASH_INT_TYPE> keys;

	TabuHasher(uint64_t seed = DEFAULT_SEED)
	{
		std::mt19937_64 rng(seed);
		for(int i=0;i<SIZEOF_HASH_INT_TYPE;i++){
			for(int j=0;j<256;j++){
				uint64_t r = rng();
				keys[i][j] = r;
			}
		}
	}


    uint64_t hash(HASH_INT_TYPE x) const
    {
		uint64_t ret = 0;
        //may likely be unlooped
        for(int i=0;i<SIZEOF_HASH_INT_TYPE;i++){
            uint8_t keyi = x & 0xff;
            ret ^= keys[i][keyi];
            x >>= 8;
        }
        return ret;
    }

    static constexpr uint64_t max()
    {
        return std::mt19937_64::max();
    }

    template <class Archive>
    void serialize(Archive &ar, const unsigned int )
    {
        ar & keys;
    }

    int64_t get_memory_usage()
    {
        return (int64_t)sizeof(uint64_t)*SIZEOF_HASH_INT_TYPE*256;
    }
};


using TabuHasher32 = TabuHasher<int32_t>;
using TabuHasher64 = TabuHasher<int64_t>;

