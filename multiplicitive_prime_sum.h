#ifndef ALGO_MULTIPLICITIVE_PRIME_SUM_H_
#define ALGO_MULTIPLICITIVE_PRIME_SUM_H_

#include <algorithm>
#include <cmath>
#include <vector>

#include "defs.h"

namespace algo {

template <typename T>
class MultiplicitivePrimeSum {
 public:
  MultiplicitivePrimeSum(llong n);

  // Calculates \sum{ f(p) : 0 < p <= M and p is a prime } for some M.
  // - mf_sum: the prefix sum function of an completely multiplicitive
  //     function f. I.e., mf_sum(n) = (\sum_{i=1}^n f(i)).
  void GetSumOverPrimes(NtFunctionT<T> mf_sum);

  // Gets the \sum{ f(p) : 0 < p <= k and p is a prime }.
  T GetSum(llong k) const {
    return k <= SG_N ? sum[k] : sum2[N/k];
  }

 private:
  void update(T& a, llong index, int pIndex) {
    int p = primes[pIndex];
    a -= (GetSum(index/p) - GetSum(p-1)) * fp[pIndex];
  }

  llong N;
  // SG_N = floor(N^{1/2})
  int SG_N;
  std::vector<int> primes;

  // Temporary varibles used during the calculation.

  // fp[i] is the f value for primes[i]. caching for better preformance
  std::vector<T> fp;

  // sum[k] = \sum{ f(p) : 0 < p <= k && p is a prime }
  // sum2[k] = \sum{ f(p) : 0 < p <= N/k && p is a prime }
  // To simplify the calculation, we define 1 to be a prime too
  std::vector<T> sum;
  std::vector<T> sum2;
};

template <typename T>
MultiplicitivePrimeSum<T>::MultiplicitivePrimeSum(llong n) : N(n) {
  SG_N = int(sqrt(n));

  auto is_prime = std::vector<bool>(SG_N+1, true);
  for (int i = 2; i <= SG_N; i++) if (is_prime[i]) {
    primes.push_back(i);
    for (llong j = 1LL*i*i; j <= SG_N; j += i) is_prime[j] = false;
  }

  fp = std::vector<T>(primes.size());
  sum = std::vector<T>(SG_N+1);
  sum2 = std::vector<T>(SG_N+1);
}

template <typename T>
void MultiplicitivePrimeSum<T>::GetSumOverPrimes(NtFunctionT<T> mf_sum) {
  for (int i = 0; i < primes.size(); i++)
    fp[i] = mf_sum(primes[i]) - mf_sum(primes[i]-1);
  for (int i = 0; i <= SG_N; i++) sum[i] = mf_sum(i);
  for (int i = 1; i <= SG_N; i++) sum2[i] = mf_sum(N/i);

  int pIndex = 0;
  for (llong p : primes) {
    // After each stage, sum_p and sum2_p contains all numbers that are:
    // - either a prime number, or
    // - a composite number whose smallest prime factor is >= primes[i]
    llong minN = p*(p-1);
    for (int j = 1; j <= SG_N && N/j > minN; j++)
      update(sum2[j], N/j, pIndex);
    for (int j = SG_N; j >= 0 && j > minN; j--)
      update(sum[j], j, pIndex);
    pIndex++;
  }
}

}  // namespace algo

#endif  // ALGO_MULTIPLICITIVE_PRIME_SUM_H_
