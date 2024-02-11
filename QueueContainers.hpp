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
    }

    /* discard contents of queue, prepare
      to start over. */
    void clear()
    {
        v.clear(); // does not change capacity
    }

    /* add values */
    template<typename IterT>
    inline
    void append(IterT begin, IterT end)
    {
        v.insert(v.end(), begin, end);
    }

    /* access values */
    inline const
    uint64_t& operator[](size_t idx) const
    {
        return v[idx];
    }

    /* container back */
    inline const
    uint64_t& back() const
    {
        return v.back();
    }

    /* container size */
    inline
    size_t size() const
    {
        return v.size();
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
    std::tuple<
        mpfr_holder,
        mpfr_holder,
        mpfr_holder,
        mpfr_holder> tmp_prods;

    size_t v_size;
    std::vector<std::tuple<
        mpfr_holder,
        mpfr_holder,
        mpfr_holder,
        mpfr_holder>> v;

    TmpProdContainer()
    {
        v.reserve(TargetPrimeQueueSize/N);
        v.resize(v.capacity());
        clear();

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

    inline
    void reset_tmp()
    {
        NumFactorsInTmp = 0;
        mpfr_set_ui(std::get<0>(tmp_prods), 1, MPFR_RNDD);
        mpfr_set_ui(std::get<2>(tmp_prods), 1, MPFR_RNDD);
    }

    /* discard contents of queue, prepare
      to start over. */
    void clear()
    {
        v_size = 0;
        reset_tmp();
    }

    /* add values */
    template<typename IterT>
    inline
    void append(IterT begin, IterT end)
    {
        while(begin != end)
        {
            auto DistToPause = std::min(
                N - NumFactorsInTmp,
                std::distance(begin, end));
            IterT PausePos = begin + DistToPause;
            NumFactorsInTmp += DistToPause;

            for(; begin != PausePos; begin++)
            {
                uint64_t& pval = *begin;
                mpfr_mul_ui_fast(
                    std::get<0>(tmp_prods),
                    pval+1,
                    MPFR_RNDD);
                mpfr_mul_ui_fast(
                    std::get<2>(tmp_prods),
                    pval,
                    MPFR_RNDD);
            }

            if(NumFactorsInTmp == N)
            {
                mpfr_mul(
                    std::get<1>(tmp_prods),
                    std::get<0>(tmp_prods),
                    ratio_ub_to_lb,
                    MPFR_RNDU);
                mpfr_mul(
                    std::get<3>(tmp_prods),
                    std::get<2>(tmp_prods),
                    ratio_ub_to_lb,
                    MPFR_RNDU);

                // Check if we failed to reserve
                // enough space.
                while(v_size >= v.size()) [[unlikely]]
                {
                    // Don't just resize, since that will
                    // reserve extra capacity.
                    v.reserve(2*v.capacity());
                    v.resize(v.capacity());
                }
            
                v[v_size] = tmp_prods;
                v_size++;
                reset_tmp();
            }
        }
    }

    /* access values */
    inline const
    std::tuple<
        mpfr_holder,
        mpfr_holder,
        mpfr_holder,
        mpfr_holder>& operator[](size_t idx) const
    {
        return v[idx];
    }

    /* container size */
    inline
    size_t size() const
    {
        return v_size;
    }
};

#endif
