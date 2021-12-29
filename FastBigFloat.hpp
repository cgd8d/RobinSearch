#ifndef FastBigFloat_hpp
#define FastBigFloat_hpp

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

    void set_mpfr_rndd(mpfr_t x)
    {
        


#endif
