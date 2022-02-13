#ifndef FastBigFloat_hpp
#define FastBigFloat_hpp

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
We guarantee sig[N-1] is non-zero unless the value is zero.
*/

static_assert(8 == sizeof(uint64_t),
    "uint64_t is not 8 bytes.");

static_assert(sizeof(unsigned long int) == sizeof(uint64_t),
    "Unsigned long int is not 64 bits.");

static_assert(sizeof(unsigned long long int) == sizeof(uint64_t),
    "Unsigned long long int is not 64 bits.");

template <size_t N>
struct FastBigFloat
{
    int64_t exp;
    std::array<uint64_t, N> sig;

    void set_ui(uint64_t x)
    {
        sig.fill(0);
        sig[N-1] = x;
        exp = 1-N;
    }

    void mul_ui_rndd(uint64_t x)
    {
        // In a majority of cases the highest word is
        // nonzero because x > 2^32.  So, we store the
        // lowest word separately and only shift it in
        // for the minority of cases when it is needed.
        uint64_t carry = 0;
        uint64_t lo = _mulx_u64(sig[0], x, reinterpret_cast<unsigned long long*>(&sig[0]));
        for(size_t i = 1; i < N; i++)
        {
            uint64_t tmp = _mulx_u64(sig[i], x, reinterpret_cast<unsigned long long*>(&sig[i]));
            sig[i-1] = __builtin_addcl(sig[i-1], tmp, carry, &carry);
        }
        sig[N-1] += carry; // Will not overflow.

        /*
        Note that now we need to conditionally move.
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
        bool DoMove = (sig[N-1] == 0);
        for(size_t i = N-1; i > 0; i--)
        {
            sig[i] = (DoMove ? sig[i-1] : sig[i]);
        }
        sig[0] = (DoMove ? lo : sig[0]);
        exp = (DoMove ? exp : exp+1);
    }

    void mul_ui_rndu(uint64_t x)
    {
        mul_ui_rndd(x);

        // Then increment by one.
        sig[0]++;
        if(__builtin_expect(sig[0] == 0, 0)) // [[unlikely]]
        {
            bool inc_exp = true;
            for(size_t i = 1; i < N; i++)
            {
                sig[i]++;
                if(sig[i] != 0)
                {
                    inc_exp = false;
                    break;
                }
            }
            if(inc_exp)
            {
                sig[N-1] = 1;
                exp++;
            }
        }
    }

    // X should already be initialized.
    void get_rndd(mpfr_t x)
    {
        mpfr_set_ui(x, sig[N-1], MPFR_RNDD);
        for(size_t i = N-1; i > 0; i--)
        {
            mpfr_mul_2si(x, x, 64, MPFR_RNDD);
            mpfr_add_ui(x, x, sig[i-1], MPFR_RNDD);
        }
        mpfr_mul_2si(x, x, 64*exp, MPFR_RNDD);
    }

    // X should already be initialized.
    void get_rndu(mpfr_t x)
    {
        mpfr_set_ui(x, sig[N-1], MPFR_RNDU);
        for(size_t i = N-1; i > 0; i--)
        {
            mpfr_mul_2si(x, x, 64, MPFR_RNDU);
            mpfr_add_ui(x, x, sig[i-1], MPFR_RNDU);
        }
        mpfr_mul_2si(x, x, 64*exp, MPFR_RNDU);
    }
};

template<size_t N>
std::ostream& operator<<(std::ostream& os, const FastBigFloat<N>& x)
{
    for(size_t i = N; i > 0; i--)
    {
        os << x.sig[i-1] << " ";
    }
    os << "x (2^64)^" << x.exp;
    return os;
}

#endif
