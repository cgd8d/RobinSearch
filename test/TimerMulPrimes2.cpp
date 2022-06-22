/// @example primesieve_iterator.cpp
/// Iterate over primes using primesieve::iterator.

#include <primesieve.hpp>
#include <mpfr.h>
#include <iostream>

const mpfr_prec_t Precision = 127;

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

int main(int argc, char** argv)
{
  uint64_t limit = 135452688167ull;

  primesieve::iterator it;
  uint64_t prime = it.next_prime();
  mpfr_t prod;
  mpfr_init2(prod, Precision);
  mpfr_set_ui(prod, 1, MPFR_RNDU);
  mpfr_t tmp;
  mpfr_init2(tmp, Precision);
  while(prime < limit)
  {
    mpfr_set_ui(tmp, prime, MPFR_RNDU); // exact
    mpfr_mul(prod, prod, tmp, MPFR_RNDD);
    prime = it.next_prime();
  }

  std::cout << "Product of primes <= " << limit << " = " << prod << std::endl;
  return 0;
}
