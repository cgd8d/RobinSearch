#include <cstdint>
#include <cstdio>
#include <stdexcept>
#include <boost/type_index.hpp>
#include <mpfr.h>

static_assert(
    boost::typeindex::type_id<mpfr_exp_t>() ==
    boost::typeindex::type_id<int64_t>(),
    "Type error: mpfr_exp_t");

const mpfr_prec_t Precision = 128;

/*
Convention: All mpfr_t values should have their rounding
mode appended, typically rndu or rndd.
*/

int main()
{
    int ret;
    mpfr_t ExpGamma_rndu;
    mpfr_init2(ExpGamma_rndu, Precision);
    mpfr_const_euler(ExpGamma_rndu, MPFR_RNDU);
    mpfr_exp(ExpGamma_rndu, ExpGamma_rndu, MPFR_RNDU);
    ret = mpfr_printf("ExpGamma <= %.20RgU\n", ExpGamma_rndu);
    if(ret < 0) throw std::runtime_error("failed output");
}
