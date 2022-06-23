/// @example primesieve_iterator.cpp
/// Iterate over primes using primesieve::iterator.

#include <primesieve.hpp>
#include <mpfr.h>
#include <iostream>
#include "../FastBigFloat.hpp"

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

  FastBigFloat<3> prod;
  prod.set_ui(1);
  mpfr_t tmp;
  mpfr_init2(tmp, Precision);
  while(prime < limit)
  {
    prod.mul_ui_rndd(prime);
    prime = it.next_prime();
  }
  prod.get_rndd(tmp);
  std::cout << "Product of primes <= " << limit << " = " << tmp << std::endl;
  return 0;
}
