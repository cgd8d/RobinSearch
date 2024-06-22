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

template<uint64_t NThreads>
uint64_t func4(uint64_t start, uint64_t stop)
{
  std::vector<uint64_t> acc(NThreads,0);
  uint64_t step_per_thread = (stop-start)/NThreads;

  #pragma omp parallel for num_threads(NThreads)
  for(int i = 0; i < NThreads; i++)
  {
    uint64_t start_local = start+i*step_per_thread;
    uint64_t stop_local = start_local+step_per_thread;
    primesieve::iterator it(start_local);
    it.generate_next_primes();

    for (; it.primes_[it.size_ - 1] < stop_local; it.generate_next_primes())
      acc[i] = std::accumulate(
        it.primes_,
        it.primes_ + it.size_,
        acc[i]);
    for (std::size_t j = 0; it.primes_[j] < stop_local; j++)
      acc[i] += it.primes_[j];
  }

  return std::accumulate(acc.begin(), acc.end(), 0ull);
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
  std::cout
    << "Prime sieve default size is "
    << primesieve::get_sieve_size()
    << " KiB (kibibytes)."
    << std::endl;

  for(uint64_t start :
    {0ull, 1ull << 46, 1ull << 50, 1ull << 54})
  {
    TimeFunc(func1, "func1", start, 1ull << 33, 1ull << 5);
    TimeFunc(func1, "func1", start, 1ull << 38, 1ull << 0);
    TimeFunc(func2, "func2", start, 1ull << 33, 1ull << 5);
    TimeFunc(func2, "func2", start, 1ull << 38, 1ull << 0);
    TimeFunc(func3, "func3", start, 1ull << 33, 1ull << 5);

    // with four-way multitasking.
    TimeFunc(func4<4>, "func4<4>", start, 1ull << 35, 1ull << 5);
    TimeFunc(func4<4>, "func4<4>", start, 1ull << 40, 1ull << 0);
    TimeFunc(func4<2>, "func4<2>", start, 1ull << 34, 1ull << 5);
    TimeFunc(func4<2>, "func4<2>", start, 1ull << 39, 1ull << 0);
  }
  return 0;
}
