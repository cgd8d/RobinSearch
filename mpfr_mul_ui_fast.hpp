#ifndef MPFR_MUL_UI_FAST_HPP
#define MPFR_MUL_UI_FAST_HPP

// Function that assumes:
// * All values are regular.
// * u >= 2
// * The size of x is known at computer time.
// * The operation is in-place.
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
... To finish later 


}

#endif // MPFR_MUL_UI_FAST_HPP
