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
        // In a majority of cases the highest word is
        // nonzero because x > 2^32.  So, we store the
        // lowest word separately and only shift it in
        // for the minority of cases when it is needed.
        uint64_t carry = 0;
        uint64_t lo = _mulx_u64(sig[0], x, &sig[0]);
        for(size_t i = 1; i < N; i++)
        {
            tmp = _mulx_u64(sig[i], x, &sig[i]);
            sig[i-1] = __builtin_addcll(tmp, sig[i-1], carry, &carry);
        }
        sig[N-1] += carry; // Will not overflow.

        // Note that now we need to conditionally move.
        // It is probably most efficient to use a branch
        // and rely on CPU branch prediction.  Overall
        // sig[N-1] will typically be zero in early parts
        // if execution because primes are small, and
        // nonzero later. Branch prediction should handle
        // this fine. By contrast, conditional moves
        // introduce a dependency (and delay) which
        // will harm the critical path.
        if(sig[N-1] == 0)
        {
            for(size_t i = N-1; i > 0; i++)
            {
                sig[










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
