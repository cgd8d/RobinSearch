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

    void mul_ui_rndd(uint64_t x)
    {
        std::array<uint64_t, N+1> tmp;
        tmp[0] = _mulx_u64(sig[0], x, &tmp[1]);
        for(size_t i = 0; i < N; i++)
        {
            unit lo 




    void set_mpfr_rndd(mpfr_t x)
    {
        size_t mpfr_numwords = (x->_mpfr_prec-1)/64+1;
        set_ui(0);
        size_t numwords_both = std::min(N, mpfr_numwords);
        for(size_t i = N - numwords_both;
            i < N;
            i++)
        {
            sig[i] = x->_mpfr_d[mpfr_numwords-numwords_both+i];
        ]
        exp = ...


#endif
