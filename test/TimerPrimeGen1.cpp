/// @example primesieve_iterator.cpp
/// Iterate over primes using primesieve::iterator.
/// This is based on libprimesieve with modifications.

#include <primesieve.hpp>
#include <iostream>
#include <vector>
#include <numeric>

template <typename T>
void store_primes_modified(uint64_t start,
                         uint64_t stop,
                         T& primes)
{
  if (start > 0)
    start--;
  if (~stop == 0)
    stop--;

  if (start < stop)
  {
    using V = typename T::value_type;
    std::size_t size = primes.size() + primesieve::prime_count_approx(start, stop);
    primes.reserve(size);

    primesieve::iterator it(start, stop);
    uint64_t prime = it.next_prime();
    for (; prime <= stop; prime = it.next_prime())
      primes.push_back((V) prime);
  }
}

int main(int argc, char** argv)
{
  uint64_t num_iter = (1 << 7);
  uint64_t step_size = (1 << 30);
  std::vector<uint64_t> vec;
  uint64_t sum = 0;

  for(uint64_t i = 0; i < num_iter; i++)
  {
    vec.clear();
    store_primes_modified(
      i*step_size,
      (i+1)*step_size,
      vec);
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
