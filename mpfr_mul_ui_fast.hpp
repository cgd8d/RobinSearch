#ifndef MPFR_MUL_UI_FAST_HPP
#define MPFR_MUL_UI_FAST_HPP

#include <mpfr.h>
#include <immintrin.h>

// Function that assumes:
// * All values are regular.
// * 2 <= u < 2^63.
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

inline
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

    // Count leading zeros.
    // out2 is guaranteed to be nonzero
    // because x is normalized and u >= 2.
    //int ls = __builtin_clzl(out2);

    // Do shift operations.
    // Note that there is a lot of pressure
    // on CPU port 1, so consider changing
    // this instruction to be implemented
    // by shl rather than shld (by changing
    // the or to a plus).  But last time
    // I checked it didn't help.
    // Note that shift by a count of 64 is
    // undefined in C++, so we also need to
    // ensure ls > 0.  This is guaranteed when
    // u < 2^63.
    //xp[0] = (out1 << ls) | (out0 >> (64-ls));
    //xp[1] = (out2 << ls) | (out1 >> (64-ls));

    // use shld.
    /*
    xp[0] = (uint64_t)(((((unsigned __int128)out1 << 64) | (unsigned __int128)out0) << (ls & 63)) >> 64);
    xp[1] = (uint64_t)(((((unsigned __int128)out2 << 64) | (unsigned __int128)out1) << (ls & 63)) >> 64);
    */
    /*
    ls &= 63;
    xp[0] = ls ? (out1 << ls) | (out0 >> (64-ls)) : out1;
    xp[1] = ls ? (out2 << ls) | (out1 >> (64-ls)) : out2;
    */

    xp[0] = out1;
    xp[1] = out2;
    asm(
        "lzcntq %[hi], %%rcx;\n"
        "shldq %%cl, %[mid_ro], %[hi];\n"
        "shldq %%cl, %[lo], %[mid];\n"
        : [hi] "+r" (xp[1]),
          [mid] "+r" (xp[0])
        : [mid_ro] "r" (out1),
          [lo] "r" (out0)
        : "rcx", "cc"
    );



        
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


/* an alternative version which interleaves
two mul operations to try to hide latency.
xa <- xa*a
xb <- xb*b
both are rounded down for simplicity.
*/
inline
void mpfr_mul_ui_fast_2way (
    mpfr_ptr xa,
    mpfr_ptr xb,
    unsigned long long int a,
    unsigned long long int b
)
{
    // Note: X data is stored as little endian

    unsigned long long int out0a, out1a, out2a;
    unsigned long long int p0a, p1a, p2a;
    unsigned char c0a;
    mp_limb_t *xpa = xa->_mpfr_d;

    unsigned long long int out0b, out1b, out2b;
    unsigned long long int p0b, p1b, p2b;
    unsigned char c0b;
    mp_limb_t *xpb = xb->_mpfr_d;

    // Do full multiplication x*u -> out.
    out0a = _mulx_u64(xpa[0], a, &p0a);
    out0b = _mulx_u64(xpb[0], b, &p0b);
    
    p1a = _mulx_u64(xpa[1], a, &p2a);
    p1b = _mulx_u64(xpb[1], b, &p2b);
    
    c0a = _addcarry_u64(0, p0a, p1a, &out1a);
    c0b = _addcarry_u64(0, p0b, p1b, &out1b);
    
    out2a = p2a + c0a;
    out2b = p2b + c0b;

    // Count leading zeros.
    // out2 is guaranteed to be nonzero
    // because x is normalized and u >= 2.
    int lsa = __builtin_clzl(out2a);
    int lsb = __builtin_clzl(out2b);

    // Do shift operations.
    // Note that there is a lot of pressure
    // on CPU port 1, so consider changing
    // this instruction to be implemented
    // by shl rather than shld (by changing
    // the or to a plus).  But last time
    // I checked it didn't help.
    // Note that shift by a count of 64 is
    // undefined in C++, so we also need to
    // ensure ls > 0.  This is guaranteed when
    // u < 2^63.
    xpa[0] = (out1a << lsa) | (out0a >> (64-lsa));
    xpa[1] = (out2a << lsa) | (out1a >> (64-lsa));
    xpb[0] = (out1b << lsb) | (out0b >> (64-lsb));
    xpb[1] = (out2b << lsb) | (out1b >> (64-lsb));

    // Update exp.
    xa->_mpfr_exp += (64-lsa);
    xb->_mpfr_exp += (64-lsb);
}

#endif // MPFR_MUL_UI_FAST_HPP
