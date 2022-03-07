#ifndef FastBigFloat_hpp
#define FastBigFloat_hpp

#include <tuple>
#include <utility>
#include <ostream>
#include <concepts>
#include <mpfr.h>
#include <boost/mp11.hpp>
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
    // Define tuple of N uint64_t
    using SingleType = std::tuple<uint64_t>;
    using TupleType = boost::mp11::mp_repeat_c<SingleType, N>;

    int64_t exp;
    TupleType sig;

    void set_ui(uint64_t x)
    {
        std::apply([](std::same_as<uint64_t> auto&... val)
            {
                (val = ... = 0);
            },
            sig);
        std::get<N-1>(sig) = x;
        exp = 1-N;
    }

    void mul_ui_rndd(uint64_t x)
    {
        // In a majority of cases the highest word is
        // nonzero because x > 2^32.  So, we store the
        // lowest word separately and only shift it in
        // for the minority of cases when it is needed.
        uint64_t carry = 0;
        uint64_t lo = _mulx_u64(std::get<0>(sig), x, reinterpret_cast<unsigned long long*>(&std::get<0>(sig)));

        // For I=0 to N-2, run helperfunc with sig[I]=lo and
        // sig[I+1]=hi.  Multiply is x*hi.
        auto helperfunc = [](uint64_t& lo, uint64_t& hi, uint64_t& carry, const uint64_t& x)
            {
                uint64_t tmp = _mulx_u64(hi, x, reinterpret_cast<unsigned long long*>(&hi));
                lo = __builtin_addcl(lo, tmp, carry, &carry);
            };
        auto helperfunc2 = [&]<std::size_t ...I>
            (std::index_sequence<I...>)
            {
                (helperfunc(std::get<I>(sig), std::get<I+1>(sig), carry, x),...);
            };
        helperfunc2(std::make_index_sequence<N-1>{});

        std::get<N-1>(sig) += carry; // Will not overflow.

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
        if(std::get<N-1>(sig) == 0)
        {
            auto helperfunc3 = []
                (uint64_t& LHS, const uint64_t& RHS)
                {
                    LHS = RHS;
                };
            auto helperfunc4 = [&]<std::size_t ...I>
                (std::index_sequence<I...>)
                {
                    (helperfunc3(std::get<N-I-1>(sig), std::get<N-I-2>(sig)),...);
                };
            helperfunc4(std::make_index_sequence<N-1>{});
            std::get<0>(sig) = lo;
        }
        else
        {
            exp++;
        }
    }

    void mul_ui_rndu(uint64_t x)
    {
        mul_ui_rndd(x);

        // Then increment by one.
        std::get<0>(sig)++;
        if(__builtin_expect(std::get<0>(sig) == 0, 0)) // [[unlikely]]
        {
            // If the LSW overflows then we need to carry.
            bool inc_exp = true;
            auto helper1 = [&]
                (uint64_t& val)
                {
                    if(inc_exp) val++;
                    if(val != 0) inc_exp = false;
                };
            std::apply([&]
                (uint64_t& first, auto&... rest)
                {
                    (helper1(rest),...);
                },
                sig);

            // This lady but is only relevant if sig was
            // 2^(64N)-1 and overflowed to the next word.
            if(inc_exp)
            {
                std::get<N-1>(sig) = 1;
                exp++;
            }
        }
    }

    // X should already be initialized.
    void get_rndd(mpfr_t x)
    {
        mpfr_set_ui(x, std::get<N-1>(sig), MPFR_RNDD);
        auto helper1 = [&]
            (const uint64_t& val)
            {
                mpfr_mul_2si(x, x, 64, MPFR_RNDD);
                mpfr_add_ui(x, x, val, MPFR_RNDD);
            };
        auto helper2 = [&]<std::size_t ...I>
            (std::index_sequence<I...>)
            {
                (helper1(std::get<N-I-2>(sig)),...);
            };
        helper2(std::make_index_sequence<N-1>{});

/*

        for(size_t i = N-1; i > 0; i--)
        {
            mpfr_mul_2si(x, x, 64, MPFR_RNDD);
            mpfr_add_ui(x, x, sig[i-1], MPFR_RNDD);
        }
*/

        mpfr_mul_2si(x, x, 64*exp, MPFR_RNDD);
    }

    // X should already be initialized.
    void get_rndu(mpfr_t x)
    {
        mpfr_set_ui(x, std::get<N-1>(sig), MPFR_RNDU);
        auto helper1 = [&]
            (const uint64_t& val)
            {
                mpfr_mul_2si(x, x, 64, MPFR_RNDU);
                mpfr_add_ui(x, x, val, MPFR_RNDU);
            };
        auto helper2 = [&]<std::size_t ...I>
            (std::index_sequence<I...>)
            {
                (helper1(std::get<N-I-2>(sig)),...);
            };
        helper2(std::make_index_sequence<N-1>{});

/*
        for(size_t i = N-1; i > 0; i--)
        {
            mpfr_mul_2si(x, x, 64, MPFR_RNDU);
            mpfr_add_ui(x, x, sig[i-1], MPFR_RNDU);
        }
*/

        mpfr_mul_2si(x, x, 64*exp, MPFR_RNDU);
    }
};

template<size_t N>
std::ostream& operator<<(std::ostream& os, const FastBigFloat<N>& x)
{
    auto helper1 = [&]
        (const uint64_t& val)
        {
            os << val << " ";
        };
        auto helper2 = [&]<std::size_t ...I>
            (std::index_sequence<I...>)
            {
                (helper1(std::get<N-I-1>(x.sig)),...);
            };
        helper2(std::make_index_sequence<N>{});

/*

    for(size_t i = N; i > 0; i--)
    {
        os << x.sig[i-1] << " ";
    }
*/


    os << "x (2^64)^" << x.exp;
    return os;
}

#endif
