/// @example primesieve_iterator.cpp
/// Iterate over primes using primesieve::iterator.
/// This is based on libprimesieve with modifications.

#include <primesieve.hpp>
#include <iostream>
#include <vector>

std::vector<uint64_t> PrimeQueue(1 << 14);

int main(int argc, char** argv)
{
  uint64_t limit = 135452688167ull;

  primesieve::iterator it_exp(0, 1000);

  primesieve::iterator it;
  uint64_t prime = it.next_prime();
  uint64_t sum = 0;
  size_t NextPrimeIdx = PrimeQueue.size();

  while(true)
  {
        size_t i = 0;
        while(i < PrimeQueue.size())
        {
            // Hack into primesieve iterator to enable
            // fast copy.
            // Handle i_ the way iterator usually does
            // so the iterator remains in a valid state.
            it.i_++;
            if(it.i_ == it.size_)
            {
                // Note generate_next_primes has post-
                // condition that i_ == 0.
                it.generate_next_primes();
            }
            size_t num_copy = std::min(
                PrimeQueue.size()-i,
                it.size_-it.i_);
            std::copy(
                it.primes_+it.i_,
                it.primes_+it.i_+num_copy,
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

          sum += PrimeQueue[i];
          i++;
       }

      // Does another prime iterator mess up the cache?
      sum += it_exp.next_prime();
    }


  return 0;
}
