#include <cstdint>
#include <cstdio>
#include <stdexcept>
#include <cassert>
#include <iostream>
#include <list>
#include <queue>
#include <stack>
#include <mpfr.h>
#include <primesieve.hpp>
#include <boost/circular_buffer.hpp>
#include "PlotDelta.hpp"
#include "FastBigFloat.hpp"

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

    static_assert(8 == sizeof(uint64_t),
        "uint64_t is not 8 bytes.");

    static_assert(sizeof(unsigned long int) == sizeof(uint64_t),
        "Unsigned long int is not 64 bits.");

    static_assert(sizeof(unsigned long long int) == sizeof(uint64_t),
        "Unsigned long long int is not 64 bits.");

    static_assert(sizeof(mp_limb_t) == sizeof(uint64_t),
        "mp_limb_t is not 64 bits.");
}

// Helper function to print mpfr_t with error checking.
std::ostream& operator<<(std::ostream& os, mpfr_t op)
{
    char* tmp_str;
    int ret;
    ret = mpfr_asprintf(
        &tmp_str,
        "%.10RDg",
        op);
    if(ret < 0)
    {
        throw std::runtime_error("mpfr_asprintf failed");
    }
    os << tmp_str;
    return os;
}

/*
Struct to produce plot of deltas.
*/
double exp_gamma;
PlotDeltaStruct PlotDelta;

/*
Counters for optimization.
*/
uint64_t cnt_NumPrimeFactors = 0;
uint64_t cnt_NumUniquePrimeFactors = 0;
uint64_t cnt_EpsEvalForExpZero = 0;
uint64_t cnt_LogLogNUpdates = 0;

/*
Struct to compute critical epsilon values and
hold temporary mpft_t values.
*/
struct ComputeEpsilonStruct
{
    mpfr_t tmp_mpfr;

    ComputeEpsilonStruct()
    {
        mpfr_init2(tmp_mpfr, Precision);
    }

    ~ComputeEpsilonStruct()
    {
        mpfr_clear(tmp_mpfr);
    }

    /*
    Critical epsilon for transition from
    p^e to p^(e+1), rounded down and up resp.
    This occurs when:
        sigma(p^e)/p^(e*(1+eps)) =
        sigma(p^(e+1))/p^((e+1)*(1+eps))
    i.e., when:
        p^eps
            = sigma(p^(e+1))/(p*sigma(p^e))
            = sigma(p^(e+1))/(sigma(p^(e+1))-1)
            = 1 + 1/(p^(e+1)+p^e+...+p)
    I think it is hard to rule out overflow if
    sigma(p^(e+1)) is computed in uint64_t, so to
    avoid any risk just do it all in mpfr_t.
    */
    void Do_rndd(mpfr_t rop, uint64_t p, uint8_t e)
    {
        mpfr_set_ui(tmp_mpfr, p, MPFR_RNDU);
        for(uint8_t i = 0; i < e; i++)
        {
            mpfr_add_ui(tmp_mpfr, tmp_mpfr, 1, MPFR_RNDU);
            mpfr_mul_ui(tmp_mpfr, tmp_mpfr, p, MPFR_RNDU);
        }
        mpfr_ui_div(rop, 1, tmp_mpfr, MPFR_RNDD);
        mpfr_log1p(rop, rop, MPFR_RNDD);
        mpfr_log_ui(tmp_mpfr, p, MPFR_RNDU);
        mpfr_div(rop, rop, tmp_mpfr, MPFR_RNDD);
    }
    void Do_rndu(mpfr_t rop, uint64_t p, uint8_t e)
    {
        mpfr_set_ui(tmp_mpfr, p, MPFR_RNDD);
        for(uint8_t i = 0; i < e; i++)
        {
            mpfr_add_ui(tmp_mpfr, tmp_mpfr, 1, MPFR_RNDD);
            mpfr_mul_ui(tmp_mpfr, tmp_mpfr, p, MPFR_RNDD);
        }
        mpfr_ui_div(rop, 1, tmp_mpfr, MPFR_RNDU);
        mpfr_log1p(rop, rop, MPFR_RNDU);
        mpfr_log_ui(tmp_mpfr, p, MPFR_RNDD);
        mpfr_div(rop, rop, tmp_mpfr, MPFR_RNDU);
    }
}
ComputeEpsilon;

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
                      << PrimeLo << "^" << int(Exp)
                      << " and "
                      << b.PrimeLo << "^" << int(b.Exp)
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
    */
    void UpdateEpsilon()
    {
        ComputeEpsilon.Do_rndd(
            CriticalEpsilon_rndd,
            PrimeLo,
            Exp);
        ComputeEpsilon.Do_rndu(
            CriticalEpsilon_rndu,
            PrimeLo,
            Exp);
    }

    PrimeGroup(const PrimeGroup  & ) = delete;
    PrimeGroup(PrimeGroup && ) = delete;
    PrimeGroup& operator=(const PrimeGroup & ) = delete;
    PrimeGroup& operator=(PrimeGroup && ) = delete;
};

std::ostream& operator<<(std::ostream& os, const PrimeGroup & pg)
{
    os << pg.PrimeLo;
    if(pg.Exp > 1)
    {
        os << "^" << int(pg.Exp);
    }
    if(pg.PrimeLo < pg.PrimeHi)
    {
        os << "..." << pg.PrimeHi;
        if(pg.Exp > 1)
        {
            os << "^" << int(pg.Exp);
        }
    }
    return os;
}

// Store the current number as an exact form (factorized)
// and mpfr interval.
std::list<PrimeGroup> Number_factors;
mpfr_t Number_rndd;
mpfr_t Number_rndu;

// Store NloglogN.  We only maintain a lower bound on this.
mpfr_t NloglogN_rndd;

// Store sigma(N)/exp(gamma).
// I just call it LHS since in my mind it is on the left.
mpfr_t LHS_rndd;
mpfr_t LHS_rndu;

// Store a priority queue of PrimeGroups, where the top
// is the next group to increment.
auto PrimeGroupItComp =
    [](typename std::list<PrimeGroup>::iterator x,
       typename std::list<PrimeGroup>::iterator y){ return *x < *y; };
std::priority_queue<
    typename std::list<PrimeGroup>::iterator,
    std::vector<typename std::list<PrimeGroup>::iterator>,
    decltype(PrimeGroupItComp)>
    PrimeGroupQueue;

/*
Assume it has already been shown that p_max should not be
increased.  Instead, a nonzero exponent should be incremented
by one.  Do that and update all relevant numbers.
*/
void IncrementExp()
{
    // Copy the top iterator, then pop it.
    std::list<PrimeGroup>::iterator top_it =
        PrimeGroupQueue.top();
    std::list<PrimeGroup>::iterator prev_it;
    PrimeGroupQueue.pop();

    // Update numbers.
    mpfr_t tmp;
    mpfr_init2(tmp, Precision);
    // For LHS, start by computing sigma(p^exp) with rounding down.
    mpfr_set_ui(tmp, 1, MPFR_RNDD);
    for(uint8_t i = 0; i < top_it->Exp; i++)
    {
        mpfr_mul_ui(tmp, tmp, top_it->PrimeLo, MPFR_RNDD);
        mpfr_add_ui(tmp, tmp, 1, MPFR_RNDD);
    }
    mpfr_div(LHS_rndu, LHS_rndu, tmp, MPFR_RNDU);
    mpfr_mul_ui(tmp, tmp, top_it->PrimeLo, MPFR_RNDD);
    mpfr_add_ui(tmp, tmp, 1, MPFR_RNDD);
    mpfr_mul(LHS_rndd, LHS_rndd, tmp, MPFR_RNDD);
    // Then compute sigma((p^exp) with rounding up.
    mpfr_set_ui(tmp, 1, MPFR_RNDU);
    for(uint8_t i = 0; i < top_it->Exp; i++)
    {
        mpfr_mul_ui(tmp, tmp, top_it->PrimeLo, MPFR_RNDU);
        mpfr_add_ui(tmp, tmp, 1, MPFR_RNDU);
    }
    mpfr_div(LHS_rndd, LHS_rndd, tmp, MPFR_RNDD);
    mpfr_mul_ui(tmp, tmp, top_it->PrimeLo, MPFR_RNDU);
    mpfr_add_ui(tmp, tmp, 1, MPFR_RNDU);
    mpfr_mul(LHS_rndu, LHS_rndu, tmp, MPFR_RNDU);
    // Also update number.
    mpfr_mul_ui(Number_rndd, Number_rndd, top_it->PrimeLo, MPFR_RNDD);
    mpfr_mul_ui(Number_rndu, Number_rndu, top_it->PrimeLo, MPFR_RNDU);
    mpfr_mul_ui(NloglogN_rndd, NloglogN_rndd, top_it->PrimeLo, MPFR_RNDD);

    // This is the least likely case, but it's just easier
    // to start by verifying whether there is a previous
    // iterator or not.
    if(top_it == Number_factors.begin())
    {
        if(top_it->PrimeLo == top_it->PrimeHi)
        {
            top_it->Exp++;
        }
        else
        {
            prev_it = Number_factors.emplace(top_it);
            prev_it->PrimeLo = prev_it->PrimeIter.next_prime();
            prev_it->PrimeHi = prev_it->PrimeLo;
            prev_it->Exp = top_it->Exp + 1;
            prev_it->UpdateEpsilon();
            top_it->PrimeLo = top_it->PrimeIter.next_prime();
            PrimeGroupQueue.push(prev_it);
        }
        top_it->UpdateEpsilon();
        PrimeGroupQueue.push(top_it);
    }
    else
    {
        prev_it = top_it;
        prev_it--;
        if(prev_it->Exp == top_it->Exp + 1)
        {
            prev_it->PrimeHi = top_it->PrimeLo;
            if(top_it->PrimeLo == top_it->PrimeHi)
            {
                Number_factors.erase(top_it);
            }
            else
            {
                top_it->PrimeLo = top_it->PrimeIter.next_prime();
                top_it->UpdateEpsilon();
                PrimeGroupQueue.push(top_it);
            }
        }
        else if(top_it->PrimeLo == top_it->PrimeHi)
        {
            top_it->Exp++;
            top_it->UpdateEpsilon();
            PrimeGroupQueue.push(top_it);
        }
        else
        {
            prev_it = Number_factors.emplace(top_it);
            prev_it->PrimeIter.skipto(top_it->PrimeLo-1);
            prev_it->PrimeLo = prev_it->PrimeIter.next_prime();
            assert(prev_it->PrimeLo == top_it->PrimeLo);
            prev_it->PrimeHi = prev_it->PrimeLo;
            prev_it->Exp = top_it->Exp + 1;
            prev_it->UpdateEpsilon();
            PrimeGroupQueue.push(prev_it);
            top_it->PrimeLo = top_it->PrimeIter.next_prime();
            top_it->UpdateEpsilon();
            PrimeGroupQueue.push(top_it);
        }
    }
}

double PrintNum_DeltaRatio = 2;
double NextPrintDelta = 1;
void CheckNumber()
{
    cnt_NumPrimeFactors++;
    if(mpfr_greaterequal_p(LHS_rndu, NloglogN_rndd))
    {
        // First, update values.
        mpfr_set(NloglogN_rndd, Number_rndd, MPFR_RNDD);
        mpfr_log(NloglogN_rndd, NloglogN_rndd, MPFR_RNDD);
        mpfr_log(NloglogN_rndd, NloglogN_rndd, MPFR_RNDD);
        double LogLogN_d = mpfr_get_d(NloglogN_rndd, MPFR_RNDD);
        mpfr_mul(NloglogN_rndd, NloglogN_rndd, Number_rndd, MPFR_RNDD);
        cnt_LogLogNUpdates++;

        // Compute delta = exp(gamma) loglogN - sigma(N)/N.
        mpfr_t tmp_mpfr;
        mpfr_init2(tmp_mpfr, Precision);
        mpfr_sub(tmp_mpfr, NloglogN_rndd, LHS_rndu, MPFR_RNDD);
        mpfr_div(tmp_mpfr, tmp_mpfr, Number_rndu, MPFR_RNDD);
        double delta = exp_gamma*mpfr_get_d(tmp_mpfr, MPFR_RNDD);
        mpfr_clear(tmp_mpfr);

        // Go ahead and print information.
        if(LogLogN_d > 2.5 and delta <= NextPrintDelta)
        {
            std::cout << "Updating logs on:" << std::endl;
            std::cout << "N = ";
            auto it = Number_factors.begin();
            while(true)
            {
                std::cout << *it;
                it++;
                if(it == Number_factors.end())
                {
                    break;
                }
                std::cout << " * ";
            }
            std::cout << std::endl;
            std::cout << "  = (" << Number_rndd
                      << ", " << Number_rndu << ")" << std::endl;
            std::cout << "exp(gamma)*loglogN - sigma(N)/N > ";
            std::cout << delta << std::endl;
            NextPrintDelta = delta/PrintNum_DeltaRatio;
        }

        if(LogLogN_d > 2.5 and delta > 0)
        {
            PlotDelta.AddPoint(LogLogN_d, std::log(delta));
        }

        // Finally, check if violation persists.
        if(mpfr_greaterequal_p(LHS_rndu, NloglogN_rndd))
        {
            std::cout << "Maybe the bound is violated?" << std::endl;
            std::cout << "NloglogN > " << NloglogN_rndd << std::endl;
            std::cout << "sigma(N)/exp(gamma) = (" << LHS_rndd << ", " << LHS_rndu << ")" << std::endl;
            if(mpfr_cmp_ui(Number_rndu, 5040) > 0)
            {
                exit(0);
            }
            else
            {
                std::cout << "But N is small so this is expected." << std::endl;
            }
        }
    }
}

/*
Add prime factors with exponent 1.
Return the number of prime factors added.
Repeat calls may add more factors until the
return value is zero.
Internally we record calculations of the critical
epsilon, since the goal is to compute it for a very
small subset of factors.
*/
boost::circular_buffer<uint64_t> PrimeQueue(4096);
primesieve::iterator PrimeQueueProducer;
struct PrimeQueueEpsilonGroup
{
    uint64_t index;
    mpfr_t Epsilon_rndu;
    PrimeQueueEpsilonGroup(uint64_t idx)
    : index(idx)
    {
        mpfr_init2(Epsilon_rndu, Precision);
        ComputeEpsilon.Do_rndu(
            Epsilon_rndu,
            PrimeQueue[idx],
            0);
    }
    ~PrimeQueueEpsilonGroup()
    {
        mpfr_clear(Epsilon_rndu);
    }
};
std::stack<PrimeQueueEpsilonGroup> PrimeQueueEpsilonStack;
uint64_t AddPrimeFactors()
{
    // Basic requirement - the last prime group should
    // have exp of one. It just makes things much simpler
    // and ought to be true.
    if(Number_factors.empty() or Number_factors.back().Exp != 1)
    {
        std::cerr << "AddPrimeFactors called with bad Number_factors." << std::endl;
        throw std::invalid_argument("AddPrimeFactors called with bad Number_factors.");
    }

    // Fill up PrimeQueue.
    assert(PrimeQueueEpsilonStack.empty() == PrimeQueue.empty());
    if(PrimeQueueEpsilonStack.empty())
    {
        while(not PrimeQueue.full())
        {
            PrimeQueue.push_front(PrimeQueueProducer.next_prime());
        }
        PrimeQueueEpsilonStack.emplace(0);
        cnt_EpsEvalForExpZero++;
    }

    // Find a safe number of primes to add.
    while(true)
    {

        // If we have already computed an epsilon
        // which shows some primes are safe to add,
        // do that.
        if(mpfr_greaterequal_p(
            PrimeQueueEpsilonStack.top().Epsilon_rndu,
            PrimeGroupQueue.top()->CriticalEpsilon_rndd))
        {
            mpfr_t eps_rndd;
            mpfr_init2(eps_rndd, Precision);
            ComputeEpsilon.Do_rndd(
                eps_rndd,
                PrimeQueue[PrimeQueueEpsilonStack.top().index],
                0);
            if(mpfr_lessequal_p(
                eps_rndd,
                PrimeGroupQueue.top()->CriticalEpsilon_rndu))
            {
                std::cerr << "Unable to compare epsilon for "
                          << PrimeQueue[PrimeQueueEpsilonStack.top().index] << "^0"
                          << " and "
                          << PrimeGroupQueue.top()->PrimeLo << "^" << int(PrimeGroupQueue.top()->Exp)
                          << std::endl;
                throw std::runtime_error("Insufficient accuracy for epsilon.");
            }
            mpfr_clear(eps_rndd);

            uint64_t this_idx = PrimeQueue.size();
            mpfr_t mpfr_temp1, mpfr_temp2;
            mpfr_init2(mpfr_temp1, Precision);
            mpfr_init2(mpfr_temp2, Precision);
            while(this_idx > PrimeQueueEpsilonStack.top().index)
            {
                if(this_idx - PrimeQueueEpsilonStack.top().index >= 16)
                {
                    FastBigFloat<3> lhs_update_rndd;
                    lhs_update_rndd.set_ui(1);
                    FastBigFloat<3> lhs_update_rndu;
                    lhs_update_rndu.set_ui(1);
                    FastBigFloat<3> rhs_update_rndd;
                    rhs_update_rndd.set_ui(1);
                    FastBigFloat<3> rhs_update_rndu;
                    rhs_update_rndu.set_ui(1);

                    for(size_t i = this_idx - 16;
                        i < this_idx;
                        i++)
                    {
                        lhs_update_rndd.mul_ui_rndd(PrimeQueue[i]+1);
                        lhs_update_rndu.mul_ui_rndu(PrimeQueue[i]+1);
                        rhs_update_rndd.mul_ui_rndd(PrimeQueue[i]);
                        rhs_update_rndu.mul_ui_rndu(PrimeQueue[i]);
                    }
                    lhs_update_rndu.get_rndu(mpfr_temp1);
                    mpfr_mul(mpfr_temp1, mpfr_temp1, LHS_rndu, MPFR_RNDU);
                    rhs_update_rndd.get_rndd(mpfr_temp2);
                    mpfr_mul(mpfr_temp2, mpfr_temp2, NloglogN_rndd, MPFR_RNDD);








                }
                this_idx--;
                uint64_t this_p = PrimeQueue[this_idx];
                mpfr_mul_ui(Number_rndd, Number_rndd, this_p, MPFR_RNDD);
                mpfr_mul_ui(Number_rndu, Number_rndu, this_p, MPFR_RNDU);
                mpfr_mul_ui(NloglogN_rndd, NloglogN_rndd, this_p, MPFR_RNDD);
                mpfr_mul_ui(LHS_rndd, LHS_rndd, this_p+1, MPFR_RNDD);
                mpfr_mul_ui(LHS_rndu, LHS_rndu, this_p+1, MPFR_RNDU);
                Number_factors.back().PrimeHi = this_p;
                cnt_NumUniquePrimeFactors++;
                CheckNumber();
            }
            mpfr_clear(mpfr_temp1);
            mpfr_clear(mpfr_temp2);

            uint64_t retval = PrimeQueue.size() - this_idx;
            PrimeQueue.erase_end(retval);
            PrimeQueueEpsilonStack.pop();
            return retval;
        }

        // Have we already shown that no primes can be
        // safely added?
        if(PrimeQueueEpsilonStack.top().index + 1 == PrimeQueue.size())
        {
            return 0;
        }

        // The last possibility is that we should compute
        // another epsilon.
        uint64_t new_idx = (PrimeQueueEpsilonStack.top().index + PrimeQueue.size())/2;
        PrimeQueueEpsilonStack.emplace(new_idx);
        cnt_EpsEvalForExpZero++;
    }
}

int main()
{
    CheckTypes();
    mpfr_init2(Number_rndd, Precision);
    mpfr_init2(Number_rndu, Precision);
    mpfr_init2(NloglogN_rndd, Precision);
    mpfr_init2(LHS_rndd, Precision);
    mpfr_init2(LHS_rndu, Precision);

    // Initialize everything to N = 1.
    mpfr_const_euler(LHS_rndd, MPFR_RNDU);
    mpfr_exp(LHS_rndd, LHS_rndd, MPFR_RNDU);
    mpfr_ui_div(LHS_rndd, 1, LHS_rndd, MPFR_RNDD);
    mpfr_const_euler(LHS_rndu, MPFR_RNDD);
    mpfr_exp(LHS_rndu, LHS_rndu, MPFR_RNDD);
    exp_gamma = mpfr_get_d(LHS_rndu, MPFR_RNDD);
    mpfr_ui_div(LHS_rndu, 1, LHS_rndu, MPFR_RNDU);
    mpfr_set_ui(Number_rndd, 1, MPFR_RNDD);
    mpfr_set_ui(Number_rndu, 1, MPFR_RNDU);
    mpfr_set_ui(NloglogN_rndd, 0, MPFR_RNDD);
    Number_factors.resize(1);
    Number_factors.front().PrimeLo =
        Number_factors.front().PrimeIter.next_prime();
    Number_factors.front().PrimeHi =
        Number_factors.front().PrimeLo;
    Number_factors.front().Exp = 0;
    Number_factors.front().UpdateEpsilon();
    PrimeGroupQueue.push(Number_factors.begin());

    // Step forward to N=2 to initiate processing.
    IncrementExp();
    PrimeQueueProducer.next_prime(); // discard value.
    cnt_NumPrimeFactors++;
    cnt_NumUniquePrimeFactors++;

    // Continue processing.
    while(Number_factors.front().Exp < 36)
    {
        while(true)
        {
            uint64_t NumFactors = AddPrimeFactors();
            if(NumFactors == 0) break;
        }
        IncrementExp();
        CheckNumber();
    }

    // Print counters.
    std::cout << "cnt_NumPrimeFactors = " << cnt_NumPrimeFactors << std::endl;
    std::cout << "cnt_NumUniquePrimeFactors = " << cnt_NumUniquePrimeFactors << std::endl;
    std::cout << "cnt_EpsEvalForExpZero = " << cnt_EpsEvalForExpZero << std::endl;
    std::cout << "cnt_LogLogNUpdates = " << cnt_LogLogNUpdates << std::endl;

    mpfr_clear(Number_rndd);
    mpfr_clear(Number_rndu);
    mpfr_clear(NloglogN_rndd);
    mpfr_clear(LHS_rndd);
    mpfr_clear(LHS_rndu);
}
