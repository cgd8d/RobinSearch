#ifndef FastBigFloat_hpp
#define FastBigFloat_hpp

/*
Support fast multiplication with multi-word precision.
The trade-off is that we don't guarantee the MSB to be
non-zero, only the most significant word.  Bad for addition,
but ok for this use case.
*/

template <size_t N>
struct FastBigFloat
{
    uint64_t exp;


#endif
