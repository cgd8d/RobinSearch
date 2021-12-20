#include <cstdint>
#include <cstdio>
#include <stdexcept>
#include <iostream>
#include <list>
#include <queue>
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

    static_assert(sizeof(unsigned long int) == sizeof(uint64_t),
        "Unsigned long int is not 64 bits.");
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
    Compare by critical epsilon value.
    This should only be called when the critical epsilon
    has been initialized.  If the two epsilons are not
    strictly ordered, raise an exception.
    x < y when x's epsilon is less than y's, meaning that
    y should have its exponent increased before x.
    */
    bool operator <(const PrimeGroup &b) const
    {
        if(mpfr_less_p(CriticalEpsilon_rndu, b.CriticalEpsilon_rndd))
        {
            return true;
        }
        else if(mpfr_greater_p(CriticalEpsilon_rndd, b.CriticalEpsilon_rndu))
        {
            return false;
        }
        else
        {
            std::cerr << "Unable to compare epsilon for "
                      << PrimeLo << "^" << Exp
                      << " and "
                      << b.PrimeLo << "^" << b.Exp
                      << std::endl;
            throw std::runtime_error("Insufficient accuracy for epsilon.");
        }
    }
    bool operator >(const PrimeGroup &b) const
    {
        return b < *this;
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
    void UpdateEpsilon()
    {
        mpfr_t temp;
        mpfr_init2(temp, Precision);

        // First compute CriticalEpsilon_rndd.
        mpfr_set_ui(temp, PrimeLo, MPFR_RNDU);
        for(uint8_t i = 0; i < Exp; i++)
        {
            mpfr_add_ui(temp, temp, 1, MPFR_RNDU);
            mpfr_mul_ui(temp, temp, 1, MPFR_RNDU);
        }
        mpfr_ui_div(CriticalEpsilon_rndd, 1, temp, MPFR_RNDD);
        mpfr_log1p(CriticalEpsilon_rndd, CriticalEpsilon_rndd, MPFR_RNDD);
        mpfr_set_ui(temp, PrimeLo, MPFR_RNDU);
        mpfr_log(temp, temp, MPFR_RNDU);
        mpfr_div(CriticalEpsilon_rndd, CriticalEpsilon_rndd, temp, MPFR_RNDD);

        // Then compute CriticalEpsilon_rndu.
        mpfr_set_ui(temp, PrimeLo, MPFR_RNDD);
        for(uint8_t i = 0; i < Exp; i++)
        {
            mpfr_add_ui(temp, temp, 1, MPFR_RNDD);
            mpfr_mul_ui(temp, temp, 1, MPFR_RNDD);
        }
        mpfr_ui_div(CriticalEpsilon_rndu, 1, temp, MPFR_RNDU);
        mpfr_log1p(CriticalEpsilon_rndu, CriticalEpsilon_rndu, MPFR_RNDU);
        mpfr_set_ui(temp, PrimeLo, MPFR_RNDD);
        mpfr_log(temp, temp, MPFR_RNDD);
        mpfr_div(CriticalEpsilon_rndu, CriticalEpsilon_rndu, temp, MPFR_RNDU);

        mpfr_clear(temp);
    }

    PrimeGroup(const PrimeGroup  & ) = delete;
    PrimeGroup(PrimeGroup && ) = delete;
    PrimeGroup& operator=(const PrimeGroup & ) = delete;
    PrimeGroup& operator=(PrimeGroup && ) = delete;
};

// Store the current number.
std::list<PrimeGroup> CurrentNumber;

// Store a priority queue of PrimeGroups, where the top
// is the next group to increment.
auto PrimeGroupItComp =
    [](typename std::list<PrimeGroup>::iterator x,
       typename std::list<PrimeGroup>::iterator y){ return *x < *y; };

std::priority_queue<
    typename std::list<PrimeGroup>::iterator,
//    std::vector<typename std::list<PrimeGroup>::iterator>,
    std::vector<auto>,
    decltype(PrimeGroupItComp)>
    PrimeGroupQueue;


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
