#ifndef ALGO_MULTIPLICITIVE_SUM_H_
#define ALGO_MULTIPLICITIVE_SUM_H_

#include <algorithm>
#include <cmath>
#include <cstring>
#include <memory>
#include <vector>

#include "defs.h"
#include "modular.h"
#include "multiplicitive_prime_sum.h"

namespace algo {

template <typename T>
class MultiplicitiveCombination {
 public:
  typedef std::tuple<llong, NtFunctionT<T>> CombinationEntry;

  class Builder {
   public:
    Builder() {}

    Builder& AddFunction(llong coef, NtFunctionT<T> func) {
      entries.push_back(std::make_tuple(coef, func));
      return *this;
    }

    MultiplicitiveCombination Build() {
      return MultiplicitiveCombination(entries);
    }
   private:
    std::vector<CombinationEntry> entries;
  };

  MultiplicitiveCombination(
      std::vector<CombinationEntry> functions_)
          : functions(functions_) {}

  std::vector<CombinationEntry> GetFunctions() const { return functions; }

 private:
  std::vector<CombinationEntry> functions;
};

template <typename T>
class MultiplicitiveSum {
 public:
  // Parameters:
  // - n: the max index to sum. i.e., all values in range [1, n].
  MultiplicitiveSum(llong n);

  // Gets the sum of given multiplicitive functions in range [1, n].
  // Parameters:
  // - mf: The multiplicitive function. mf(p, e) = f(p^e),
  //     where p will always be a prime when calling the function.
  //     f(1) should always be 1, otherwise all f(n) = 0.
  // - mf_sum: the prefix sum function of an approximate fully
  //     multiplicitive function g. I.e.,
  //     mf_sum(n) = (\sum_{i=1}^n g(i)).
  //     Here g must be a fully multiplicitive function, and
  //     g(p) = f(p) for all prime p.
  T PrefixSum(MultiplicitiveFunctionT<T> mf,
              NtFunctionT<T> mf_sum);

  // Same prefix function, but uses linear combination of fully
  // multiplicitive functions to approximate the targt function.
  T PrefixSum(MultiplicitiveFunctionT<T> mf,
              const MultiplicitiveCombination<T>& mc);

 private:
  // Calculates \sum{ f(p) : 0 < p <= M and p is a prime } for some M.
  void GetSumOverPrimes(const MultiplicitiveCombination<T>& mc);

  // Calculates the sum of all f(m) where:
  // - m is a composite number, and
  // - m is a multiple of X (not including given X), and
  // - X < m < N
  // Parameters:
  // - X: the given number to find multiples
  // - fx: the function value f(X)
  // - pIndex: The index of max prime factors of X (-1 if X = 1)
  T GetSumOverMultiples(MultiplicitiveFunctionT<T> mf,
                        llong X, T fx, int pIndex) const;

  T getSumP(llong k) const {
    return k <= SG_N ? sum[k] : sum2[N/k];
  }

  llong N;
  // SG_N = floor(N^{1/2})
  int SG_N;
  std::vector<int> primes;

  // Temporary varibles used during the calculation.

  // sum[k] = \sum{ f(p) : 0 < p <= k && p is a prime }
  // sum2[k] = \sum{ f(p) : 0 < p <= N/k && p is a prime }
  // To simplify the calculation, we define 1 to be a prime too
  std::vector<T> sum;
  std::vector<T> sum2;
};

template <typename T>
MultiplicitiveSum<T>::MultiplicitiveSum(llong n) : N(n) {
  SG_N = int(sqrt(n));

  auto is_prime = std::vector<bool>(SG_N+1, true);

  for (int i = 2; i <= SG_N; i++) if (is_prime[i]) {
    primes.push_back(i);
    for (llong j = 1LL*i*i; j <= SG_N; j += i) is_prime[j] = false;
  }
}

template <typename T>
void MultiplicitiveSum<T>::GetSumOverPrimes(
    const MultiplicitiveCombination<T>& mc) {
  sum = std::vector<T>(SG_N+1, T(0));
  sum2 = std::vector<T>(SG_N+1, T(0));

  MultiplicitivePrimeSum<T> mps(N);
  for (auto entry : mc.GetFunctions()) {
    mps.GetSumOverPrimes(std::get<1>(entry));
    llong coef = std::get<0>(entry);
    for (int i = 1; i <= SG_N; i++) {
      sum[i] += coef*mps.GetSum(i);
      sum2[i] += coef*mps.GetSum(N/i);
    }
  }
}

template <typename T>
T MultiplicitiveSum<T>::GetSumOverMultiples(
    MultiplicitiveFunctionT<T> mf, llong X, T fx, int pIndex) const {
  T res = 0;

  // Categorize the multiples Y into two buckets and calculate them separately:
  // 1. Y has only 1 largest prime factor
  // 2. Y has >= 2 largest prime factors
  for (int i = pIndex+1; i < primes.size(); i++) {
    if (X == N) printf("%d/%d\n", i, primes.size());
    int p = primes[i];
    if (1LL*p*p > X) break;

    llong nextX = X;
    for (int j = 1; nextX >= p; j++) {
      nextX /= p;
      T nextF = fx * mf(p, j);
      if (nextX > p) {
        // Summarize numbers in 1st category
        llong sumP = getSumP(nextX) - getSumP(p);
        res += nextF * sumP;
        res += GetSumOverMultiples(mf, nextX, nextF, i);
      }
      // Summarize numbers in 2nd category
      if (j >= 2) res += nextF;
    }
  }
  return res;
}

template <typename T>
T MultiplicitiveSum<T>::PrefixSum(MultiplicitiveFunctionT<T> mf,
                                  NtFunctionT<T> mf_sum) {
  typename MultiplicitiveCombination<T>::Builder builder;

  builder.AddFunction(1, mf_sum);
  return PrefixSum(mf, builder.Build());
}

template <typename T>
T MultiplicitiveSum<T>::PrefixSum(MultiplicitiveFunctionT<T> mf,
                                  const MultiplicitiveCombination<T>& mc) {
  GetSumOverPrimes(mc);

  T res = GetSumOverMultiples(mf, N, T(1), -1) + getSumP(N) - getSumP(1) + 1;
  return res;
}

}  // namespace algo

#endif  // ALGO_MULTIPLICITIVE_SUM_H_
