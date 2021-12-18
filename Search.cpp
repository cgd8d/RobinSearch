#include <cstdint>
#include <cstdio>
#include <stdexcept>
#include <limits>
#include <mpfr.h>

const mpfr_prec_t Precision = 128;

/*
Convention: All mpfr_t values should have their rounding
mode appended, typically rndu or rndd.
*/

void CheckTypes()
{
    if(mpfr_get_emax() != std::numeric_limits<int64_t>::max())
    {
        throw std::logic_error("mpfr_get_emax is unexpected");
    }
}


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
