#ifndef QueueContainers_hpp
#define QueueContainers_hpp

/*
The purpose of this file is to hold queued
data precomputed for use in adding new primes.
This is in two parts:
* PrimeQueueContainer, to hold primes that have
not been added yet.
* TmpProdContainer, to hold products of bunches
of values.

These containers are special-purpose because
they are essential for performance and have
access patterns which I may be able to take
advantage of.  Putting it into a structure
makes it easier to experiment with that.
*/

#include <vector>
#include <tuple>

const size_t TargetPrimeQueueSize = 1 << 26;

struct PrimeQueueContainer
{

    std::vector<uint64_t> v;

    PrimeQueueContainer()
    {
        v.reserve(TargetPrimeQueueSize);
    };

    /* discard contents of queue, prepare
      to start over. */
    void reset()
    {
        v.clear(); // does not change capacity
    }

    /* add values */
    inline
    template<typename IterT>
    void append(IterT begin, IterT end)
    {
        v.insert(v.end(), begin, end)
    }

    /* access values */
    inline
    uint64_t operator[](size_t idx) const
    {
        return v[idx];
    }

};

template<uint64_t N> 
struct TmpProdContainer
{



    std::vector<std::tuple<
        mpfr_holder,
        mpfr_holder,
        mpfr_holder,
        mpfr_holder>> v;

    TmpProdContainer()
    {
        v.reserve(TargetPrimeQueueSize/N);
    }


        
};

#endif
