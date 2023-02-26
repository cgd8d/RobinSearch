#ifndef MPFR_MUL_UI_FAST_HPP
#define MPFR_MUL_UI_FAST_HPP

#include <iostream>
#include <mpfr.h>
#include <immintrin.h>

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

void mpfr_mul_ui_fast (mpfr_ptr x, unsigned long long int u, mpfr_rnd_t rnd_mode)
{
    // Note: X data is stored as little endian

    unsigned long long int out0, out1, out2;
    unsigned long long int p0, p1, p2;
    unsigned char c0;
    mp_limb_t *xp = x->_mpfr_d;

    // Do full multiplication x*u -> out.
    out0 = _mulx_u64(xp[0], u, &p0);
    p1 = _mulx_u64(xp[1], u, &p2);
    c0 = _addcarry_u64(0, p0, p1, &out1);
    out2 = p2 + c0;

    std::cout << "(" << xp[1] << ", " << xp[0] << ") * " << u << std::endl;
    std::cout << " = (" << out2 << ", " << out1 << ", " << out0 << ")" << std::endl;

    // Count leading zeros.
    // out2 is guaranteed to be nonzero
    // because x is normalized and u >= 2.
    int ls = __builtin_clzl(out2);
    std::cout << "ls = " << ls << std::endl;

    // Do shift operations.
    xp[0] = (out1 << ls) | (out0 >> (64-ls));
    xp[1] = (out2 << ls) | (out1 >> (64-ls));

    std::cout << "out2 << " << ls << " = " << (out2<<ls) << std::endl;
    std::cout << "after shift, x is (" << xp[1] << ", " << xp[0] << ") << std::endl;


    // Update exp.
    x->_mpfr_exp += (64-ls);

    // Rounding.
    // Slightly conservative since the discarded
    // bits could be zero, but we accept that
    // for the time savings of not checking.
    if(rnd_mode == MPFR_RNDU)
    {
        if(xp[0] == (unsigned long int)(-1)) [[unlikely]]
        {
            mpfr_nextabove(x);
        }
        else
        {
            xp[0]++;
        }
    }
}

#endif // MPFR_MUL_UI_FAST_HPP
