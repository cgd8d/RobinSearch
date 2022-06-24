/// @example primesieve_iterator.cpp
/// Iterate over primes using primesieve::iterator.
/// This is based on libprimesieve with modifications.

#include <primesieve.hpp>
#include <iostream>
#include <vector>
#include <numeric>
#include <algorithm>

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
    while(true)
    {
      it.i_++;
      if(it.i_ >= it.size_)
        it.generate_next_primes();

      uint64_t* it_begin = it.primes_ + it.i_;
      uint64_t* it_last = it.primes_ + it.size_;

/*
      uint64_t* it_end;
      if(it.primes_[it.size_-1] < stop)
        it_end = it_last;
      else
        it_end = std::lower_bound(it_begin, it_last, stop);
      primes.insert(primes.end(), it_begin, it_end);
      it.i_ += std::distance(it_begin, it_end)-1;
      if(it_end < it_last)
        break;
*/
      if(it.primes_[it.size_-1] < stop)
      {
        primes.insert(primes.end(), it_begin, it_end);
        it.i_ = it_last-1;
      } else {

      


    }
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
