/// @example primesieve_iterator.cpp
/// Iterate over primes using primesieve::iterator.
/// This is based on libprimesieve with modifications.

#include <primesieve.hpp>
#include <mpfr.h>
#include <iostream>
#include <vector>

const mpfr_prec_t Precision = 128;

std::vector<uint64_t> PrimeQueue(1 << 14);

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

  primesieve::iterator it_exp(0, 1000);

  primesieve::iterator it;
  uint64_t prime = it.next_prime();
  mpfr_t sum;
  mpfr_init2(sum, Precision);
  mpfr_set_ui(sum, 1, MPFR_RNDU);
  size_t NextPrimeIdx = PrimeQueue.size();

  mpfr_t tmp;
  mpfr_init2(tmp, Precision);

  while(true)
  {
        size_t i = 0;
        while(i < PrimeQueue.size())
        {
            // Hack into primesieve iterator to enable
            // fast copy.
            // Handle i_ the way iterator usually does
            // so the iterator remains in a valid state.
            if(it.i_++ == it.last_idx_)
            {
                // Note generate_next_primes has post-
                // condition that i_ == 0.
                it.generate_next_primes();
            }
            size_t num_copy = std::min(
                PrimeQueue.size()-i,
                it.last_idx_+1-it.i_);
            std::copy(
                it.primes_.begin()+it.i_,
                it.primes_.begin()+it.i_+num_copy,
                PrimeQueue.begin()+i);
            it.i_ += num_copy-1;
            i += num_copy;
        }
        NextPrimeIdx = 0;


        i = 0;
        while(i < PrimeQueue.size())
        {
          if(PrimeQueue[i] > limit)
          {
            std::cout << "Sum of primes <= " << limit << " = " << sum << std::endl;

            // Note that since sum is a 64-bit variable the result
            // will be incorrect (due to integer overflow) if
            // limit > 10^10. However we do allow limits > 10^10
            // since this is useful for benchmarking.
            if (limit > 10000000000ull)
              std::cerr << "Warning: sum is likely incorrect due to 64-bit integer overflow!" << std::endl;

            return 0;
          }

          mpfr_set_ui(tmp, PrimeQueue[i], MPFR_RNDU); // exact
          mpfr_mul(sum, sum, tmp, MPFR_RNDD);
          i++;
       }

      // Does another prime iterator mess up the cache?
      mpfr_set_ui(tmp, it_exp.next_prime(), MPFR_RNDU); // exact
      mpfr_mul(sum, sum, tmp, MPFR_RNDD);
    }


  return 0;
}
