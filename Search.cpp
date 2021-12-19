#include <cstdint>
#include <cstdio>
#include <stdexcept>
#include <iostream>
#include <mpfr.h>
#include <primesieve.hpp>

const mpfr_prec_t Precision = 128;

/*
Convention: All mpfr_t values should have their rounding
mode appended, typically rndu or rndd.
*/

void CheckTypes()
{
    mpfr_set_emax(mpfr_get_emax_max());
    if(mpfr_get_emax() != (1LL << 62) - 1)
    {
        std::cerr << "mpfr_get_emax() = " << mpfr_get_emax() << std::endl;
        std::cerr << "mpfr_get_emax_max() = " << mpfr_get_emax_max() << std::endl;
        throw std::logic_error("mpfr_get_emax is unexpected");
    }
}

/*
Structure to store a group of prime factors with
the same exponent.  The range PrimeLo to PrimeHi
is inclusive.
PrimeIter.next_prime() should always return the
prime just after PrimeLo.
*/
struct PrimeGroup
{
    uint64_t PrimeLo;
    uint64_t PrimeHi;
    uint8_t Exp;
    primesieve::iterator PrimeIter;
    mpfr_t CriticalEpsilon_rndd;
    mpfr_t CriticalEpsilon_rndu;

    /*
    Initialize the PrimeIter with a modest stop_hint,
    which cuts down on memory usage for the majority
    of instantiations which will not iterate very high.
    */
    PrimeGroup()
    : PrimeIter(0, 1000)
    {
        mpfr_init2(CriticalEpsilon_rndd, Precision);
        mpfr_init2(CriticalEpsilon_rndu, Precision);
    }

    ~PrimeGroup()
    {
        mpfr_clear(CriticalEpsilon_rndd);
        mpfr_clear(CriticalEpsilon_rndu);
    }

    /*
    Compute the critical epsilon value that leads to
    incrementing the exponent of PrimeLo from Exp to Exp+1.
    This occurs when:
        sigma(PrimeLo^Exp)/PrimeLo^(Exp*(1+eps)) =
        sigma(PrimeLo^(Exp+1))/PrimeLo^((Exp+1)*(1+eps))
    i.e., when:
        PrimeLo^eps
            = sigma(PrimeLo^(Exp+1))/(PrimeLo*sigma(PrimeLo^Exp))
            = sigma(PrimeLo^(Exp+1))/(sigma(PrimeLo^(Exp+1))-1)
            = 1 + 1/(PrimeLo^(Exp+1)+PrimeLo^Exp+...+PrimeLo)
    I think it is hard to rule out overflow if
    sigma(PrimeLo^(Exp+1)) is computed in uint64_t, so to
    avoid any risk just do it all in mpfr_t.
    */
};


int main()
{
    CheckTypes();
    int ret;
    mpfr_t ExpGamma_rndu;
    mpfr_init2(ExpGamma_rndu, Precision);
    mpfr_const_euler(ExpGamma_rndu, MPFR_RNDU);
    mpfr_exp(ExpGamma_rndu, ExpGamma_rndu, MPFR_RNDU);
    ret = mpfr_printf("ExpGamma <= %.20RgU\n", ExpGamma_rndu);
    if(ret < 0) throw std::runtime_error("failed output");
}
