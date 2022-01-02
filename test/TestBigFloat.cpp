#include <random>
#include <upstream>
#include <exception>
#include <mpfr.h>
#include "../BigFloat.hpp"

void CheckIntervals(BigFloat<3>& bf1, BigFloat<3>& bf2,
                    mpfr_t mp1, mpfr_t mp2)
{
    mpfr_t t1;
    mpfr_init2(t1, 1024);
    bf1.get_rndd(t1);
    if(mpfr_greater_p(t1, mp2))
    {
        throw std::runtime_error("Intervals do not overlap.");
    }
    mpfr_t t2;
    mpfr_init2(t2, 1024);
    bf2.get_rndu(t2);
    if(mpfr_less_p(t2, mp1))
    {
        throw std::runtime_error("Intervals do not overlap.");
    }
    mpfr_sub(t1, t2, t1, MPFR_RNDD);
    mpfr_sub(t2, mp2, mp1, MPFR_RNDD);
    mpfr_div(t1, t1, t2, MPFR_RNDD);
    if(mpfr_less_p(mp1, mp2) and mpfr_cmp_ui(t1, 10) > 0)
    {
        throw std::runtime_error("BigFloat interval is too big.");
    }
}

int main()
{

    std::mt19937_64 mt;

    // Test starting at 1.
    for(size_t i = 0; i < 1000; i++)
    {
        BigFloat<3> bf1, bf2;
        mpfr_t mp1, mp2;
        bf1.set_ui(1);
        bf2.set_ui(1);
        mpfr_init2(mp1, 128);
        mpfr_init2(mp2, 128);
        mpfr_set_ui(mp1, 1, MPFR_RNDD);
        mpfr_set_ui(mp2, 1, MPFR_RNDU);
        CheckIntervals(bf1, bf2, mp1, mp2);
        for(size_t j = 0; j < 1000; j++)
        {
            uint64_t x = mt();
            bf1.mul_ui_rndd(x);
            bf2.mul_ui_rndu(x);
            mpfr_mul_ui(mp1, mp1, x, MPFR_RNDD);
            mpfr_mul_ui(mp2, mp2, x, MPFR_RNDU);
            CheckIntervals(bf1, bf2, mp1, mp2);
        }
    }

// Next, should test starting from arbitrary.
    std::cout << "passed test of BigFloat." << std::endl;
}
