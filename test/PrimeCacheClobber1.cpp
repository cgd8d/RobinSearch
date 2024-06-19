/// @example primesieve_iterator.cpp
/// Iterate over primes using primesieve::iterator.
/// This is based on libprimesieve with modifications.

#include <primesieve.hpp>
#include <iostream>
#include <chrono>

uint64_t func1(uint64_t start, uint64_t stop)
{
  uint64_t acc = 0
  primesieve::iterator it(start);
  uint64_t prime = it.next_prime();
  for (; prime < stop; prime = it.next_prime())
    acc += prime;
  return acc;
}

template<typename F>
void TimeFunc(F func, std::string fname)
{
  uint64_t step_per_chunk = (1ull << 27);
  uint64_t nchunks = (1ull << 10);
  uint64_t sum_mod64 = 0;
  
  auto t1 = std::high_resolution_clock::now();
  for(size_t i = 0; i < nchunks; i++) {
    uint64_t start = i*step_per_chunk;
    uint64_t stop = start+step_per_chunk;
    sum_mod64 += func(start, stop);
  }
  auto t2 = std::high_resolution_clock::now();

  std::cout << "Using " << fname
    << ", sum of primes < " << step_per_chunk*nchunks
    << " (mod 2^64) = " << sum_mod64
    << " computed in " << t2-t1
    << std::endl;
}

int main(int argc, char** argv)
{
  TimeFunc(func1, "func1");

  return 0;
}
