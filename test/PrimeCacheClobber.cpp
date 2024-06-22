/// @example primesieve_iterator.cpp
/// Iterate over primes using primesieve::iterator.
/// This is based on libprimesieve with modifications.

#include <primesieve.hpp>
#include <iostream>
#include <chrono>
#include <numeric>
using namespace std::chrono_literals;

uint64_t func1(uint64_t start, uint64_t stop)
{
  uint64_t acc = 0;
  primesieve::iterator it(start);
  it.generate_next_primes();

  for (; it.primes_[it.size_ - 1] < stop; it.generate_next_primes())
    acc = std::accumulate(
      it.primes_,
      it.primes_ + it.size_,
      acc);
  for (std::size_t i = 0; it.primes_[i] < stop; i++)
    acc += it.primes_[i];

  return acc;
}

uint64_t func2(uint64_t start, uint64_t stop)
{
  uint64_t acc = 0;
  primesieve::iterator it(start);
  uint64_t prime = it.next_prime();
  for (; prime < stop; prime = it.next_prime())
    acc += prime;
  return acc;
}

uint64_t func3(uint64_t start, uint64_t stop)
{
  std::vector<uint64_t> primes;
  primesieve::generate_primes(start, stop, &primes);
  uint64_t acc = std::accumulate(
    primes.begin(),
    primes.end(),
    0ull);
  return acc;
}

template<typename F>
void TimeFunc(
    F func,
    std::string fname,
    uint64_t start,
    uint64_t step_per_chunk,
    uint64_t nchunks
)
{
  uint64_t sum_mod64 = 0;
  
  auto t1 = std::chrono::high_resolution_clock::now();
  for(size_t i = 0; i < nchunks; i++) {
    uint64_t this_start = start+i*step_per_chunk;
    uint64_t this_stop = this_start+step_per_chunk;
    sum_mod64 += func(this_start, this_stop);
  }
  auto t2 = std::chrono::high_resolution_clock::now();

  std::cout << "Using " << fname
    << ", sum of primes in ["
    << start << ", " << start+step_per_chunk*nchunks
    << "] (mod 2^64) = " << sum_mod64
    << " computed with " << nchunks
    << " chunks in " << (t2-t1)/1.0s << " s"
    << std::endl;
}

/*
Test assuming prime queue length of 2^28
which roughly corresponds to a step size
of 2^33.
*/

int main(int argc, char** argv)
{
  for(uint64_t start :
    {0ull, 1ull << 46, 1ull << 50, 1ull << 54})
  {
    TimeFunc(func1, "func1", start, 1ull << 33, 1ull << 5);
    TimeFunc(func1, "func1", start, 1ull << 38, 1ull << 0);
    TimeFunc(func2, "func2", start, 1ull << 33, 1ull << 5);
    TimeFunc(func2, "func2", start, 1ull << 38, 1ull << 0);
    TimeFunc(func3, "func3", start, 1ull << 33, 1ull << 5);
  }
  return 0;
}
