#ifndef Robin_mpfr_helpers_hpp
#define Robin_mpfr_helpers_hpp

#include <mpfr.h>

const mpfr_prec_t Precision = 128;
const size_t NumLimbs = 2;

// Helper object to store initialized mpfr objects.
// Provide correct copy semantics.
// Also provide implicit conversion to mpfr_t.
struct mpfr_holder
{
    MPFR_DECL_INIT(val, Precision);

    mpfr_holder() = default;
    mpfr_holder(const mpfr_holder& src)
    {
        val->_mpfr_prec = Precision;
        val->_mpfr_sign = 1;
        //*val = *src.val;
        val->_mpfr_d = __gmpfr_local_tab_val;
        *this = src;
        /*for(size_t i = 0;
            i < NumLimbs;
            i++)
        {
            __gmpfr_local_tab_val[i] = src.__gmpfr_local_tab_val[i];
        }*/
    }
    // assume this is already a well-formed object.
    // also assume it is positive.
    mpfr_holder& operator=(const mpfr_holder& src)
    {
        //*val = *src.val;
        //val->_mpfr_d = __gmpfr_local_tab_val;
        val->_mpfr_exp = src.val->_mpfr_exp;
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

#endif
