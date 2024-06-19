/// @example primesieve_iterator.cpp
/// Iterate over primes using primesieve::iterator.
/// This is based on libprimesieve with modifications.

#include <primesieve.hpp>
#include <iostream>
#include <vector>
#include <numeric>

int main(int argc, char** argv)
{
  uint64_t step_per_chunk = (1ull << 27);
  uint64_t nchunks = (1ull << 10);
  uint64_t sum_mod64 = 0;

  for(size_t i = 0; i < nchunks; i++) {
    uint64_t start = i*step_per_chunk;
    uint64_t end = start+step_per_chunk;
    std::vector<uint64_t> primes;
    primesieve::generate_primes(start, end, &primes);
    sum_mod64 += std::accumulate(
      primes.begin(),
      primes.end(),
      sum_mod64);
  }
  
  std::cout << "Sum of primes < " << step_per_chunk*nchunks << " (mod 2^64) = " << sum_mod64 << std::endl;

  return 0;
}
