/// @example primesieve_iterator.cpp
/// Iterate over primes using primesieve::iterator.
/// This is based on libprimesieve with modifications.

#include <primesieve.hpp>
#include <iostream>

int main(int argc, char** argv)
{
  uint64_t step_per_chunk = (1ull << 27);
  uint64_t nchunks = (1ull << 10);
  uint64_t sum_mod64 = 0;

  for(size_t i = 0; i < nchunks; i++) {
    uint64_t start = i*step_per_chunk;
    uint64_t end = start+step_per_chunk;
    primesieve::iterator it(start);

    uint64_t prime = it.next_prime();
    for (; prime < limit; prime = it.next_prime())
      sum_mod64 += prime;
  }
  
  std::cout << "Sum of primes < " << step_per_chunk*nchunks << " (mod 2^64) = " << sum_mod64 << std::endl;

  return 0;
}
