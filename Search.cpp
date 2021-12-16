#include <typeinfo>
#include "mpfr.h"

static_assert(typeid(mpfr_exp_t) == typeid(int64_t),
    "Type error: mpfr_exp_t");
