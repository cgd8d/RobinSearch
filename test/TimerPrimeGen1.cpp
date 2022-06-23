/// @example primesieve_iterator.cpp
/// Iterate over primes using primesieve::iterator.
/// This is based on libprimesieve with modifications.

#include <primesieve.hpp>
#include <iostream>
#include <vector>
#include <algorithm>

int main(int argc, char** argv)
{
  uint64_t num_iter = (1 << 7);
  uint64_t step_size = (1 << 30);
  std::vector<uint64_t> vec;
  primesieve::iterator it;
  uint64_t sum = 0;

  for(uint64_t i = 0; i < num_iter; i++)
  {
    vec.clear();
    primesieve::generate_primes(
      i*step_size,
      (i+1)*step_size,
      &vec);
    sum = std::accumulate(vec.begin(), vec.end(), sum);
  }

  std::cout << "Sum of primes <= " << num_iter*step_size << " = " << sum << std::endl;

  // Note that since sum is a 64-bit variable the result
  // will be incorrect (due to integer overflow) if
  // limit > 10^10. However we do allow limits > 10^10
  // since this is useful for benchmarking.
  if (num_iter*step_size > 10000000000ull)
    std::cerr << "Warning: sum is likely incorrect due to 64-bit integer overflow!" << std::endl;

  return 0;
}
