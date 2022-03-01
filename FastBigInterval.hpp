#ifndef FastBigInterval_hpp
#define FastBigInterval_hpp

#include <array>
#include <ostream>
#include <mpfr.h>
#include <immintrin.h>

/*
Support fast multiplication with multi-word precision.
The trade-off is that we don't guarantee the MSB to be
non-zero, only the most significant word.  Bad for addition,
but ok for this use case.
Also this struct represents strictly non-negative values.
The struct represents:
(sig[0] + sig[1]*2^64 + ...)*2^(64*exp)
The interval is (sig_lo, sig_hi).
We guarantee sig_hi[N-1] is non-zero unless the value is zero.
*/

static_assert(8 == sizeof(uint64_t),
    "uint64_t is not 8 bytes.");

static_assert(sizeof(unsigned long int) == sizeof(uint64_t),
    "Unsigned long int is not 64 bits.");

static_assert(sizeof(unsigned long long int) == sizeof(uint64_t),
    "Unsigned long long int is not 64 bits.");

template <size_t N>
struct FastBigInterval
{
    int64_t exp;
    std::array<uint64_t, N> sig_lo;
    std::array<uint64_t, N> sig_hi;

    void set_ui(uint64_t x)
    {
        sig_lo.fill(0);
        sig_lo[N-1] = x;
        sig_hi = sig_lo;
        exp = 1-N;
    }

    void mul_ui(uint64_t x)
    {
        // In a majority of cases the highest word is
        // nonzero because x > 2^32.  So, we store the
        // lowest word separately and only shift it in
        // for the minority of cases when it is needed.

        // Start with the hi value, since its msw determines
        // whether we need to shift.
        uint64_t carry_hi = 0;
        uint64_t lsw_hi = _mulx_u64(sig_hi[0], x, reinterpret_cast<unsigned long long*>(&sig_hi[0]));
        for(size_t i = 1; i < N; i++)
        {
            uint64_t tmp = _mulx_u64(sig_hi[i], x, reinterpret_cast<unsigned long long*>(&sig_hi[i]));
            sig_hi[i-1] = __builtin_addcl(sig_hi[i-1], tmp, carry_hi, &carry_hi);
        }
        sig_hi[N-1] += carry_hi; // Will not overflow.

        // Follow with lo value while we wait for hi compute
        // to settle.
        uint64_t carry_lo = 0;
        uint64_t lsw_lo = _mulx_u64(sig_lo[0], x, reinterpret_cast<unsigned long long*>(&sig_lo[0]));
        for(size_t i = 1; i < N; i++)
        {
            uint64_t tmp = _mulx_u64(sig_lo[i], x, reinterpret_cast<unsigned long long*>(&sig_lo[i]));
            sig_lo[i-1] = __builtin_addcl(sig_lo[i-1], tmp, carry_lo, &carry_lo);
        }
        sig_lo[N-1] += carry_lo; // Will not overflow.

        /*
        Now we need to conditionally move.
        When x is sufficiently large, using a branch and
        letting branch prediction work is probably best,
        but when x occupies not much more than 32 bits
        then the branch prediction won't be very reliable
        and a conditional move will probably be better.
        (We can hide the latency by running multiple
        multiplies together, since we need four per
        iteration.)
        Plausibly we could put the words into a vector
        register and shuffle, conditionally.
        Another option would be to manage sig in a circular
        buffer with capacity N+1, and rather than shifting
        words we can shift the start pointer.
        Anyway, come back and experiment with this one.
        */
        if(sig_hi[N-1] == 0)
        {
            for(size_t i = N-1; i > 0; i--)
            {
                sig_hi[i] = sig_hi[i-1];
                sig_lo[i] = sig_lo[i-1];
            }
            sig_hi[0] = lsw_hi;
            sig_lo[0] = lsw_lo;
        }
        else
        {
            assert(sig_lo[N-1] == 0);
            exp++;

            // The calculation is not exact.
            // By default sig_lo is rounded down.
            // We need to explicitly round sig_hi up.
            sig_hi[0]++;

            if(sig_hi[0] == 0) [[unlikely]]
            {
                bool inc_exp = true;
                for(size_t i = 1; i < N; i++)
                {
                    sig_hi[i]++;
                    if(sig_hi[i] != 0)
                    {
                        inc_exp = false;
                        break;
                    }
                }
                if(inc_exp)
                {
                    sig_hi[N-1] = 1;
                    exp++;
                }
            }
        }
    }

    // x_lo and x_hi should already be initialized.
    void get(mpfr_t x_lo, mpfr_t x_hi)
    {
        mpfr_set_ui(x_lo, sig_lo[N-1], MPFR_RNDD);
        for(size_t i = N-1; i > 0; i--)
        {
            mpfr_mul_2si(x_lo, x_lo, 64, MPFR_RNDD);
            mpfr_add_ui(x_lo, x_lo, sig_lo[i-1], MPFR_RNDD);
        }
        mpfr_mul_2si(x_lo, x_lo, 64*exp, MPFR_RNDD);

        mpfr_set_ui(x_hi, sig_hi[N-1], MPFR_RNDU);
        for(size_t i = N-1; i > 0; i--)
        {
            mpfr_mul_2si(x_hi, x_hi, 64, MPFR_RNDU);
            mpfr_add_ui(x_hi, x_hi, sig_hi[i-1], MPFR_RNDU);
        }
        mpfr_mul_2si(x_hi, x_hi, 64*exp, MPFR_RNDU);
    }
};

template<size_t N>
std::ostream& operator<<(std::ostream& os, const FastBigInterval<N>& x)
{
    os << "(";
    for(size_t i = N; i > 0; i--)
    {
        os << x.sig_lo[i-1] << " ";
    }
    os << "- ";
    for(size_t i = N; i > 0; i--)
    {
        os << x.sig_hi[i-1] << " ";
    }
    os << ") x (2^64)^" << x.exp;
    return os;
}


#endif
