#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <stdexcept>
#include <cassert>
#include <ctime>
#include <iostream>
#include <list>
#include <queue>
#include <stack>
#include <array>
#include <tuple>
#include <fstream>
#include <mpfr.h>
#include <omp.h>
#include <primesieve.hpp>
#include <primesieve/CpuInfo.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/list.hpp>
#include "PlotDelta.hpp"
#include "mpfr_mul_ui_fast.hpp"

const mpfr_prec_t Precision = 128;
const size_t NumLimbs = 2;
const int NumThreads = 2;

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

    if(NumThreads != omp_get_max_threads())
    {
        std::cerr << "NumThreads = " << NumThreads << std::endl;
        std::cerr << "omp_get_max_threads() = " << omp_get_max_threads() << std::endl;
        throw std::logic_error("NumThreads is wrong.");
    }

    static_assert(8 == sizeof(uint64_t),
        "uint64_t is not 8 bytes.");

    static_assert(sizeof(unsigned long int) == sizeof(uint64_t),
        "Unsigned long int is not 64 bits.");

    static_assert(sizeof(unsigned long long int) == sizeof(uint64_t),
        "Unsigned long long int is not 64 bits.");

    static_assert(sizeof(mp_limb_t) == sizeof(uint64_t),
        "mp_limb_t is not 64 bits.");

    static_assert(Precision == 64*NumLimbs,
        "Precision and NumLimbs are not consistent.");
}

// Helper object to store initialized mpfr objects.
// Provide correct copy semantics.
// Also provide implicit conversion to mpfr_t.
struct mpfr_holder
{
    MPFR_DECL_INIT(val, Precision);

    mpfr_holder() = default;
    mpfr_holder(const mpfr_holder& src)
    {
        *val = *src.val;
        val->_mpfr_d = __gmpfr_local_tab_val;
        for(size_t i = 0;
            i < NumLimbs;
            i++)
        {
            __gmpfr_local_tab_val[i] = src.__gmpfr_local_tab_val[i];
        }
    }
    mpfr_holder& operator=(mpfr_holder& src)
    {
        *val = *src.val;
        val->_mpfr_d = __gmpfr_local_tab_val;
        for(size_t i = 0;
            i < NumLimbs;
            i++)
        {
            __gmpfr_local_tab_val[i] = src.__gmpfr_local_tab_val[i];
        }
        return *this;
    }

    operator mpfr_t&()
    {
        return val;
    }

    operator const mpfr_t&() const
    {
        return val;
    }

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        // I'm not completely sure what the default
        // constructor catches, so being conservative here.
        ar & val->_mpfr_prec; // might be redundant
        ar & val->_mpfr_exp;
        ar & val->_mpfr_sign;
        ar & __gmpfr_local_tab_val;
        val->_mpfr_d = __gmpfr_local_tab_val; // might be redundant...
    }
};
std::vector<mpfr_holder> mpfr_tmp;

// Helper function to print mpfr_t with error checking.
std::ostream& operator<<(std::ostream& os, mpfr_t op)
{
    char* tmp_str;
    int ret;
    ret = mpfr_asprintf(
        &tmp_str,
        "%.40RNg",
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
uint64_t cnt_FastBunchMul = 0;
uint64_t cnt_SlowMulExpOne = 0;

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
    void Do_rndd(mpfr_t& rop, uint64_t p, uint8_t e)
    {
        mpfr_set_ui(tmp_mpfr, p, MPFR_RNDU);
        for(uint8_t i = 0; i < e; i++)
        {
            mpfr_add_ui(tmp_mpfr, tmp_mpfr, 1, MPFR_RNDU);
            mpfr_mul_ui(tmp_mpfr, tmp_mpfr, p, MPFR_RNDU);
        }
        mpfr_ui_div(rop, 1, tmp_mpfr, MPFR_RNDD);
        mpfr_log1p(rop, rop, MPFR_RNDD);
        mpfr_set_ui(tmp_mpfr, p, MPFR_RNDU);
        mpfr_log(tmp_mpfr, tmp_mpfr, MPFR_RNDU);
        mpfr_div(rop, rop, tmp_mpfr, MPFR_RNDD);
    }
    void Do_rndu(mpfr_t& rop, uint64_t p, uint8_t e)
    {
        mpfr_set_ui(tmp_mpfr, p, MPFR_RNDD);
        for(uint8_t i = 0; i < e; i++)
        {
            mpfr_add_ui(tmp_mpfr, tmp_mpfr, 1, MPFR_RNDD);
            mpfr_mul_ui(tmp_mpfr, tmp_mpfr, p, MPFR_RNDD);
        }
        mpfr_ui_div(rop, 1, tmp_mpfr, MPFR_RNDU);
        mpfr_log1p(rop, rop, MPFR_RNDU);
        mpfr_set_ui(tmp_mpfr, p, MPFR_RNDD);
        mpfr_log(tmp_mpfr, tmp_mpfr, MPFR_RNDD);
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
    
    // mutable: hackish way to check current
    // status so we can save it.
    // better would be to track last retrieved value.
    // to revisit.
    primesieve::iterator PrimeIter;
    mpfr_holder CriticalEpsilon_rndd;
    mpfr_holder CriticalEpsilon_rndu;

    /*
    Initialize the PrimeIter with a modest stop_hint,
    which cuts down on memory usage for the majority
    of instantiations which will not iterate very high.
    */
    PrimeGroup()
    : PrimeIter(0, 1000)
    {


    }

    ~PrimeGroup()
    {


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
    
    template<class Archive>
    void save(Archive & ar, const unsigned int version) const
    {
        // note, version is always the latest when saving
        ar & PrimeLo;
        ar & PrimeHi;
        ar & Exp;
        ar & CriticalEpsilon_rndd;
        ar & CriticalEpsilon_rndu;
    }
    template<class Archive>
    void load(Archive & ar, const unsigned int version)
    {
        ar & PrimeLo;
        ar & PrimeHi;
        ar & Exp;
        ar & CriticalEpsilon_rndd;
        ar & CriticalEpsilon_rndu;
        if(PrimeLo < 500)
        {
            // provide hint that this will probably start small.
            PrimeIter.jump_to(PrimeLo+1, 1000);
        }
        else
        {
            PrimeIter.jump_to(PrimeLo+1);
        }
    }
    BOOST_SERIALIZATION_SPLIT_MEMBER()
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
mpfr_holder Number_rndd;
mpfr_holder Number_rndu;

// Store NloglogN.  We only maintain a lower bound on this.
mpfr_holder NloglogN_rndd;

// Store sigma(N)/exp(gamma).
// I just call it LHS since in my mind it is on the left.
mpfr_holder LHS_rndd;
mpfr_holder LHS_rndu;

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
            prev_it->PrimeIter.jump_to(top_it->PrimeLo);
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

// Print the current number without any computation
// or checks.
void PrintNumber()
{
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
}

// Deal with checking the values and logging as appropriate.
// Return true when the logarithms are recomputed,
// false otherwise; this helps other parts of the code
// know when certain work is appropriate to redo.
double PrintNum_DeltaRatio = 2;
double NextPrintDelta = 1;
bool CheckNumber()
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
        // Acquire mpfr_tmp[0].
        mpfr_sub(mpfr_tmp[0].val, NloglogN_rndd, LHS_rndu, MPFR_RNDD);
        mpfr_div(mpfr_tmp[0].val, mpfr_tmp[0].val, Number_rndu, MPFR_RNDD);
        double delta = exp_gamma*mpfr_get_d(mpfr_tmp[0].val, MPFR_RNDD);
        // Release mpfr_tmp[0].

        // Go ahead and print information.
        if(LogLogN_d > 2.5 and delta <= NextPrintDelta)
        {
            std::cout << "Updating logs on:" << std::endl;
            PrintNumber();
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
            if(mpfr_cmp_d(Number_rndu, 5040.5) > 0)
            {
                std::cout << "And it seems to be the real deal." << std::endl;
                std::cout << "Number_rndd = " << mpfr_get_d(Number_rndd, MPFR_RNDD) << std::endl;
                std::cout << "Number_rndu = " << mpfr_get_d(Number_rndu, MPFR_RNDU) << std::endl;
                exit(0);
            }
            else
            {
                std::cout << "But N is small so this is expected." << std::endl;
            }
        }
        return true;
    }
    else
    {
        // No update to logarithms, so return false.
        return false;
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
Large queues are needed because successive calls
to primesieve::generate_primes have a startup
cost proportional to sqrt of the prime values.
Since primes will easily exceed 1e12, a sqrt(p)
penalty is a lot.
Azure free GitHub runners have 7GB RAM, see
https://docs.github.com/en/actions/using-github-hosted-runners/about-github-hosted-runners#supported-runners-and-hardware-resources
*/
std::vector<std::vector<uint64_t>> PrimeQueueVec(NumThreads);
const size_t TargetPrimeQueueSize = 1 << 28;
size_t PrimeQueueStep = 2*TargetPrimeQueueSize;
uint64_t NextPrimeToGen = 3;
size_t PrimeQueueVecIdx = NumThreads-1;
size_t NextPrimeIdx = 0;
struct PrimeQueueEpsilonGroup
{
    uint64_t index;
    mpfr_holder Epsilon_rndu;

    PrimeQueueEpsilonGroup(uint64_t idx)
    : index(idx)
    {
        ComputeEpsilon.Do_rndu(
            Epsilon_rndu,
            PrimeQueueVec[PrimeQueueVecIdx][idx],
            0);
    }
};
std::stack<PrimeQueueEpsilonGroup> PrimeQueueEpsilonStack;

// Store intermediate products of groups of primes.
const size_t ProductGroupSize = 256;
std::array<std::vector<std::tuple<
    mpfr_holder,
    mpfr_holder,
    mpfr_holder,
    mpfr_holder>>, NumThreads> TmpProducts;

// For groups of values to multiply, rather than
// separately multiplying with rounding up
// and down, multiply the whole group with
// rounding down and then compute a conservative
// upper bound with a scaling factor.
// This approach roughly doubles interval size
// but gives significant speedup.
mpfr_holder ratio_ub_to_lb;

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

    // Basic requirement - check that prime queue
    // and epsilon stack are in sync.
    assert(PrimeQueueEpsilonStack.empty() ==
           (NextPrimeIdx == PrimeQueueVec[PrimeQueueVecIdx].size()));

    // Set up PrimeQueue.
    if(NextPrimeIdx == PrimeQueueVec[PrimeQueueVecIdx].size())
    {
        // Reached the end of this prime queue,
        // so advance one.
        PrimeQueueVecIdx++;
        NextPrimeIdx = 0;
    }
    if(PrimeQueueVecIdx >= NumThreads)
    {
        // Need to generate more primes.
        #pragma omp parallel for num_threads(NumThreads)
        for(int i = 0; i < NumThreads; i++)
        {
            PrimeQueueVec[i].clear();
            PrimeQueueVec[i].reserve(TargetPrimeQueueSize);
            primesieve::generate_primes(
                NextPrimeToGen + i*PrimeQueueStep,
                NextPrimeToGen + (i+1)*PrimeQueueStep - 1,
                &PrimeQueueVec[i]);

            // Within threads, also compute
            // intermediate products.
            TmpProducts[i].clear();
            TmpProducts[i].reserve(
                TargetPrimeQueueSize/ProductGroupSize);
            TmpProducts[i].resize(
                PrimeQueueVec[i].size()/ProductGroupSize);
            for(size_t j = 0;
                j < TmpProducts[i].size();
                j++)
            {
                // Initialize lhs.
                mpfr_set_ui(
                    std::get<0>(TmpProducts[i][j]),
                    PrimeQueueVec[i][j*ProductGroupSize]+1,
                    MPFR_RNDD);

                // Initialize rhs.
                mpfr_set_ui(
                    std::get<2>(TmpProducts[i][j]),
                    PrimeQueueVec[i][j*ProductGroupSize],
                    MPFR_RNDD);

                // Iterate on multiplication.
                for(size_t k = 1;
                    k < ProductGroupSize;
                    k++)
                {
                    uint64_t pval = PrimeQueueVec[i][j*ProductGroupSize+k];
                    mpfr_mul_ui_fast(
                        std::get<0>(TmpProducts[i][j]),
                        pval+1,
                        MPFR_RNDD);
                    mpfr_mul_ui_fast(
                        std::get<2>(TmpProducts[i][j]),
                        pval,
                        MPFR_RNDD);
                }
                // Compute upper bounds based
                // on lower bounds.
                mpfr_mul(
                    std::get<1>(TmpProducts[i][j]),
                    std::get<0>(TmpProducts[i][j]),
                    ratio_ub_to_lb,
                    MPFR_RNDU);
                mpfr_mul(
                    std::get<3>(TmpProducts[i][j]),
                    std::get<2>(TmpProducts[i][j]),
                    ratio_ub_to_lb,
                    MPFR_RNDU);
            } // End loop over group j.
        } // End parallel region.
        if(PrimeQueueVec.back().back() > (1ull<<63))
        {
            std::cout << "we can't use prime "
                "values more than 2^63 because "
                "of the implementation of "
                "mpfr_mul_ui_fast." << std::endl;
            exit(0);
        }
        NextPrimeToGen += NumThreads*PrimeQueueStep;
        if(2*PrimeQueueVec.back().size() <
           TargetPrimeQueueSize)
        {
            // Next time use a bigger step.
            PrimeQueueStep *= 2;
        }
        PrimeQueueVecIdx = 0;
        assert(NextPrimeIdx == 0);
    }
    std::vector<uint64_t>& PrimeQueue = PrimeQueueVec[PrimeQueueVecIdx];
    std::vector<std::tuple<
        mpfr_holder,
        mpfr_holder,
        mpfr_holder,
        mpfr_holder>>& ThisTempProd = TmpProducts[PrimeQueueVecIdx];
    if(PrimeQueueEpsilonStack.empty())
    {
        PrimeQueueEpsilonStack.emplace(PrimeQueue.size() - 1);
        cnt_EpsEvalForExpZero++;
    }

    // Find a safe number of primes to add.
    // We want to ensure the set of primes to
    // add is maximal, meaning either that the
    // next action will be to increase an exponent
    // more than one or that we hit the end of
    // a prime queue.
    // EndPrimeToAdd is one past the last
    // prime to add.
    size_t EndPrimeToAdd = NextPrimeIdx;
    while((not PrimeQueueEpsilonStack.empty()) and
          mpfr_greaterequal_p(
              PrimeQueueEpsilonStack.top().Epsilon_rndu,
              PrimeGroupQueue.top()->CriticalEpsilon_rndd))
    {
        // PrimeQueueEpsilonStack.top() is safe
        // to add, but it's not necessarily
        // the last one safe to add.
        EndPrimeToAdd = PrimeQueueEpsilonStack.top().index+1;
        PrimeQueueEpsilonStack.pop();
    }
    // Now we have a lower bound (possibly still
    // NextPrimeIdx which would add no primes)
    // and either an upper bound from
    // PrimeQueueEpsilonStack.top() or we've
    // popped the whole stack.
    // If the stack is empty then we're done.
    if(PrimeQueueEpsilonStack.empty())
    {
        assert(EndPrimeToAdd == PrimeQueue.size());
    }
    else
    {
        size_t EndPrimeToAdd_ub = PrimeQueueEpsilonStack.top().index;
        while(EndPrimeToAdd < EndPrimeToAdd_ub)
        {
            size_t Trial_idx = (EndPrimeToAdd+EndPrimeToAdd_ub)/2;
            PrimeQueueEpsilonGroup Trial_eps(Trial_idx);
            cnt_EpsEvalForExpZero++;
            if(mpfr_greaterequal_p(
                Trial_eps.Epsilon_rndu,
                PrimeGroupQueue.top()->CriticalEpsilon_rndd))
            {
                // Trial_idx is safe to add.
                EndPrimeToAdd = Trial_idx+1;
            }
            else
            {
                // Trial_idx is definitely not safe to add.
                PrimeQueueEpsilonStack.push(Trial_eps);
                EndPrimeToAdd_ub = Trial_idx;
            }
        }
    }

    // Detect unlikely case where we can't add anything.
    if(NextPrimeIdx == EndPrimeToAdd)
    {
        return 0;
    }

    // Double-check that we are computing epsilon
    // with enough accuracy to be sure about
    // ordering.
    // Acquire mpfr_tmp[0] (holds eps_rndd)
    ComputeEpsilon.Do_rndd(
        mpfr_tmp[0].val,
        PrimeQueue[EndPrimeToAdd-1],
        0);
    if(mpfr_lessequal_p(
        mpfr_tmp[0].val,
        PrimeGroupQueue.top()->CriticalEpsilon_rndu))
    {
        std::cerr << "Unable to compare epsilon for "
                  << PrimeQueue[PrimeQueueEpsilonStack.top().index] << "^0"
                  << " and "
                  << PrimeGroupQueue.top()->PrimeLo << "^" << int(PrimeGroupQueue.top()->Exp)
                  << std::endl;
        throw std::runtime_error("Insufficient accuracy for epsilon.");
    }
    // Release mpfr_tmp[0]

    // Iterate and do multiplication.
    // Acquire mpfr_tmp
    uint64_t NextPrimeIdx_init = NextPrimeIdx;
    bool ForceSingleMul = false;
    while(NextPrimeIdx < EndPrimeToAdd)
    {
        // Iterate

        // Iterate factor by factor until we get
        // to use precomputed bunches.
        while(NextPrimeIdx < EndPrimeToAdd and
            (NextPrimeIdx%ProductGroupSize != 0 or
             ForceSingleMul or
             NextPrimeIdx + ProductGroupSize > EndPrimeToAdd))
        {
            uint64_t this_p = PrimeQueue[NextPrimeIdx];
            mpfr_mul_ui_fast(Number_rndd, this_p, MPFR_RNDD);
            mpfr_mul_ui_fast(Number_rndu, this_p, MPFR_RNDU);
            mpfr_mul_ui_fast(NloglogN_rndd, this_p, MPFR_RNDD);
            mpfr_mul_ui_fast(LHS_rndd, this_p+1, MPFR_RNDD);
            mpfr_mul_ui_fast(LHS_rndu, this_p+1, MPFR_RNDU);
            Number_factors.back().PrimeHi = this_p;
            NextPrimeIdx++;
            cnt_NumUniquePrimeFactors++;
            CheckNumber();
            // We don't care if logs were
            // recomputed.
            ForceSingleMul = false;
            cnt_SlowMulExpOne++;
        }

        // If we stopped because we're done, break.
        if(NextPrimeIdx == EndPrimeToAdd)
        {
            break;
        }

        // Use precomputed bunches.
        // Stop when we are going past the end
        // or when we need to update logs.
        // We do get one shot at updating logs
        // on bunch boundaries, but that isn't
        // guaranteed to work.
        bool LogsAreUpdated = false;
        while(NextPrimeIdx + ProductGroupSize <= EndPrimeToAdd)
        {
            // Check if test values indicate possible violation of bound.
            // Compute updated lhs rounded up.
            mpfr_mul(
                mpfr_tmp[4].val,
                std::get<1>(ThisTempProd[NextPrimeIdx/ProductGroupSize]),
                LHS_rndu,
                MPFR_RNDU);
            // Compute updated rhs rounded down.
            mpfr_mul(
                mpfr_tmp[5].val,
                std::get<2>(ThisTempProd[NextPrimeIdx/ProductGroupSize]),
                NloglogN_rndd,
                MPFR_RNDD);

            if(mpfr_less_p(mpfr_tmp[4].val, mpfr_tmp[5].val))
            {
                // LHS < RHS is guaranteed.
                // Save current progress and keep going.
                mpfr_mul(LHS_rndd,
                         LHS_rndd,
                         std::get<0>(ThisTempProd[NextPrimeIdx/ProductGroupSize]),
                         MPFR_RNDD);
                mpfr_set(LHS_rndu,
                         mpfr_tmp[4].val,
                         MPFR_RNDU);
                mpfr_mul(Number_rndd,
                         Number_rndd,
                         std::get<2>(ThisTempProd[NextPrimeIdx/ProductGroupSize]),
                         MPFR_RNDD);
                mpfr_mul(Number_rndu,
                         Number_rndu,
                         std::get<3>(ThisTempProd[NextPrimeIdx/ProductGroupSize]),
                         MPFR_RNDU);
                mpfr_set(NloglogN_rndd,
                         mpfr_tmp[5].val,
                         MPFR_RNDD);
                NextPrimeIdx += ProductGroupSize;
                cnt_NumUniquePrimeFactors += ProductGroupSize;
                Number_factors.back().PrimeHi = PrimeQueue[NextPrimeIdx-1];
                cnt_NumPrimeFactors += ProductGroupSize;
                cnt_FastBunchMul += ProductGroupSize;

                LogsAreUpdated = false;
            }
            else if(not LogsAreUpdated)
            {
                // Update logs and try one more time.
                // The update is handled by Check number.
                CheckNumber();
                cnt_NumPrimeFactors--;
                LogsAreUpdated = true;
            }
            else
            {
                // Unable to proceed so fall
                // back to incremental muls.
                ForceSingleMul = true;
                break;
            }
        } // End loop over bunches.
    } // End iterating muls.

    uint64_t retval = NextPrimeIdx - NextPrimeIdx_init;
    return retval;
}

// serialize all necessary entries
// in this program.  Depending on the
// type of archive provided, this could
// be reading from or writing to file.
// some of these aren't strictly needed
// but help size processing correctly
// from the beginning (like PrimeQueueStep).
template<class Archive>
void DoSerializeAll(Archive& ar)
{
    ar & PlotDelta;
    ar & cnt_NumPrimeFactors;
    ar & cnt_NumUniquePrimeFactors;
    ar & cnt_EpsEvalForExpZero;
    ar & cnt_LogLogNUpdates;
    ar & cnt_FastBunchMul;
    ar & cnt_SlowMulExpOne;
    ar & Number_factors;
    ar & Number_rndd;
    ar & Number_rndu;
    ar & NloglogN_rndd; // not strictly necessary
    ar & LHS_rndd;
    ar & LHS_rndu;
    ar & PrintNum_DeltaRatio;
    ar & NextPrintDelta;
    ar & PrimeQueueStep;
}

// argv[1] is max exp.
// argv[2] is start time in sec.
// argv[3] is in filename for results.
// argv[4] is out filename for results.
int main(int argc, char *argv[])
{
    if(argc < 2 or argc > 5)
    {
        std::cerr << "Incorrect number of command-line arguments: "
                  << argc
                  << std::endl;
        return 1;
    }
    unsigned long MaxExp = std::strtoul(argv[1], nullptr, 0);
    if(MaxExp >= 64)
    {
        std::cerr << "MaxExp is too large ("
                  << MaxExp
                  << ")"
                  << std::endl;
        return 1;
    }
    unsigned long end_time = -1;
    if(argc >= 3)
    {
        end_time = std::strtoul(argv[2], nullptr, 0);
        std::cout << "Job end time will be "
                  << end_time
                  << std::endl;
    }

    // Print library versions.
    std::cout << "Primesieve version is "
              << primesieve::primesieve_version()
              << std::endl;
    std::cout << "mpfr version is "
              << mpfr_get_version()
              << std::endl;

    // Sieve size check.
    std::cout << "Prime sieve default size is "
              << primesieve::get_sieve_size()
              << " KiB (kibibytes)."
              << std::endl;
    std::cout << "Primesieve reports cache sizes:"
              << std::endl
              << primesieve::cpuInfo.l1CacheBytes()
              << " bytes (L1d)"
              << std::endl
              << primesieve::cpuInfo.l2CacheBytes()
              << " bytes (L2)"
              << std::endl
              << primesieve::cpuInfo.l3CacheBytes()
              << " bytes (L3)"
              << std::endl;
    primesieve::set_sieve_size(1<<10);
    std::cout << "Prime sieve size changed to "
              << primesieve::get_sieve_size()
              << " KiB (kibibytes)."
              << std::endl;

    // Number of threads.
    std::cout << "Max number of threads: "
              << NumThreads
              << std::endl;

    CheckTypes();
    mpfr_tmp.resize(6+4*NumThreads);

    // Compute the ratio between upper and
    // lower bound for a product of
    // ProductGroupSize values.
    {
        mpfr_holder tmp_1pluseps;
        mpfr_set_ui(ratio_ub_to_lb, 1, MPFR_RNDU);
        mpfr_set_ui(tmp_1pluseps, 1, MPFR_RNDU);
        mpfr_nextabove(tmp_1pluseps);
        for(size_t i= 0; i < ProductGroupSize; i++)
        {
            mpfr_mul(
                ratio_ub_to_lb,
                ratio_ub_to_lb,
                tmp_1pluseps,
                MPFR_RNDU);
        }
    }

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
    cnt_NumPrimeFactors++;
    cnt_NumUniquePrimeFactors++;

    if(argc >= 4 and
        not std::string(argv[3]).empty())
    {
        // read values from file.
        // this will overwrite all of the
        // initiations we just did,
        // as intended.
        std::ifstream ifs(argv[3]);
        boost::archive::text_iarchive ia(ifs);
        DoSerializeAll(ia);

        // We also need to manually set up
        // the queue for new primes.
        NextPrimeToGen = Number_factors.back().PrimeHi+1;

        // We already need to regenerate the
        // prime group priority queue. serialization
        // cannot handle iterators.
        while(not PrimeGroupQueue.empty())
        {
            PrimeGroupQueue.pop();
        }
        for(auto it = Number_factors.begin();
            it != Number_factors.end();
            it++)
        {
            PrimeGroupQueue.push(it);
        }

        std::cout << "Resuming work with status:" << std::endl;
        PrintNumber();
    }

    // Continue processing.
    while(Number_factors.front().Exp < MaxExp)
    {
        while(true)
        {
            uint64_t NumFactors = AddPrimeFactors();
            if(time(nullptr) > end_time)
            {
                std::cout << "out of time at "
                          << time(nullptr)
                          << ", exiting."
                          << std::endl;
                CheckNumber();
                goto exit_loop;
            }
            if(NumFactors == 0) break;
        }
        IncrementExp();
        CheckNumber();
    }
    exit_loop:

    std::cout << "Ending process with status:" << std::endl;
    PrintNumber();

    // Print counters.
    std::cout << "cnt_NumPrimeFactors = " << cnt_NumPrimeFactors << std::endl;
    std::cout << "cnt_NumUniquePrimeFactors = " << cnt_NumUniquePrimeFactors << std::endl;
    std::cout << "cnt_EpsEvalForExpZero = " << cnt_EpsEvalForExpZero << std::endl;
    std::cout << "cnt_LogLogNUpdates = " << cnt_LogLogNUpdates << std::endl;
    std::cout << "cnt_FastBunchMul = " << cnt_FastBunchMul << std::endl;
    std::cout << "cnt_SlowMulExpOne = " << cnt_SlowMulExpOne << std::endl;

    // save info to file if requested.
    if(argc >= 5 and
        not std::string(argv[4]).empty())
    {
        std::cout << "Saving progress to: "
            << argv[4] << std::endl;
        std::ofstream ofs(argv[4]);
        boost::archive::text_oarchive oa(ofs);
        DoSerializeAll(oa);
    }
}
