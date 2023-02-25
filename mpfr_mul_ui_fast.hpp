#ifndef MPFR_MUL_UI_FAST_HPP
#define MPFR_MUL_UI_FAST_HPP

// Function that assumes:
// * All values are regular.
// * u >= 2
// * The size of x is known at computer time.
// * The operation is in-place.
// * Rounding is either up or down.
// This allows us to skip some checks in the
// standard implementation if this function.

// Verify types.
static_assert(8 == sizeof(uint64_t),
    "uint64_t is not 8 bytes.");
static_assert(sizeof(unsigned long int) == sizeof(uint64_t),
    "Unsigned long int is not 64 bits.");
static_assert(sizeof(unsigned long long int) == sizeof(uint64_t),
    "Unsigned long long int is not 64 bits.");

// tuple enables better register usage than array
// For small compile time sizes.
// Use a "magic" type for tuple<t,t,...,t>.
#include <tuple>
template <size_t N, typename Head, typename... T>
struct magic {
    using tuple_type = typename magic<N - 1, Head, Head, T...>::tuple_type;
};
template <typename... T>
struct magic<1, T...> {
    using tuple_type = std::tuple<T...>;
};


template <size_t N>
int mpfr_mul_ui_fast (mpfr_ptr x, unsigned long int u, mpfr_rnd_t rnd_mode)
{
    magic<N+1, unsigned long int>::tuple_type{} tmp;

    // X data is stored as little endian

    unsigned long int out0, out1, out2;
    unsigned long int p0, p1, p2;
    unsigned char c0;
    mp_limb_t *xp = MPFR_MANT(x);

    // Do full multiplication x*u -> out.
    out0 = _mulx_u64(xp[0], u, &p0);
    p1 = _mulx_u64(xp[1], u, &p2);
    c0 = _addcarry_u64(0, p0, p1, &out1);
    out2 = p2 + c0;

    // Count leading zeros.
    // out2 is guaranteed to be nonzero
    // because x is normalized and u >= 2.
    int ls = __builtin_clzl(out2);

    // Do shift operations.
    xp[0] = (out1 << ls) | (out0 >> (64-ls));
    xp[1] = (out2 << ls) | (out1 >> (64-ls));

    // Rounding.
    if(rnd_mode == MPFR_RNDU)
    {
        if(xp[0] == unsigned long int(-1)) [[unlikely]]

... To finish later 

Rounding, including very unlikely carry.
Adjust exp.

}

#endif // MPFR_MUL_UI_FAST_HPP
