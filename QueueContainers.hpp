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
#include "Robin_mpfr_helpers.hpp"

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
    // track how many factors have been added
    // to tmp.
    uint64_t NumFactorsInTmp;

    // For groups of values to multiply, rather than
    // separately multiplying with rounding up
    // and down, multiply the whole group with
    // rounding down and then compute a conservative
    // upper bound with a scaling factor.
    // This approach roughly doubles interval size
    // but gives significant speedup.
    mpfr_holder ratio_ub_to_lb;

    // Local storage for work on product.
    mpfr_holder tmp_lhs_lb;
    mpfr_holder tmp_rhs_lb;

    std::vector<std::tuple<
        mpfr_holder,
        mpfr_holder,
        mpfr_holder,
        mpfr_holder>> v;

    TmpProdContainer()
    {
        v.reserve(TargetPrimeQueueSize/N);
        NumFactorsInTmp = 0;
        mpfr_set_ui(tmp_lhs_lb, 1, MPFR_RNDD);
        mpfr_set_ui(tmp_rhs_lb, 1, MPFR_RNDD);

        // Compute the ratio between upper and
        // lower bound for a product of
        // ProductGroupSize values.
        mpfr_holder tmp_1pluseps;
        mpfr_set_ui(ratio_ub_to_lb, 1, MPFR_RNDU);
        mpfr_set_ui(tmp_1pluseps, 1, MPFR_RNDU);
        mpfr_nextabove(tmp_1pluseps);
        for(size_t i = 0; i < N; i++)
        {
            mpfr_mul(
                ratio_ub_to_lb,
                ratio_ub_to_lb,
                tmp_1pluseps,
                MPFR_RNDU);
        }
    }


        
};

#endif
