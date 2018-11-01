#ifndef ALGO_MULTIPLICITIVE_SUM_H_
#define ALGO_MULTIPLICITIVE_SUM_H_

#include <algorithm>
#include <functional>
#include <memory>
#include <vector>

#include "defs.h"
#include "modular.h"

namespace algo {

typedef std::function<llong(int, int, int)> MultiplicitiveFunction;
typedef std::function<llong(llong, int)> SumFunction;

class MultiplicitiveCombination {
 public:
  typedef std::function<llong(llong, int)> SumFunction;
  typedef std::tuple<llong, SumFunction> CombinationEntry;

  class Builder {
   public:
    Builder() {}

    Builder* AddFunction(llong coef, SumFunction func) {
      entries.push_back(std::make_tuple(coef, func));
      return this;
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

class MultiplicitiveSum {
 public:
  // Parameters:
  // - n: the max index to sum. i.e., all values in range [1, n].
  MultiplicitiveSum(llong n);

  // Gets the sum of given multiplicitive functions in range [1, n].
  // Parameters:
  // - mf: The multiplicitive function. mf(p, e, mod) = f(p^e) % mod,
  //     where p will always be a prime when calling the function.
  //     f(1) should always be 1, otherwise all f(n) = 0.
  // - mf_sum: the prefix sum function of an approximate fully
  //     multiplicitive function g. I.e.,
  //     mf_sum(n, mod) = (\sum_{i=1}^n g(i)) % mod.
  //     Here g must be a fully multiplicitive function, and
  //     g(p) = f(p) for all prime p.
  // - mod: 0 if no mod is required (i.e., the real value).
  llong PrefixSum(MultiplicitiveFunction mf,
                  SumFunction mf_sum, int mod);

  // Same prefix function, but uses linear combination of fully
  // multiplicitive functions to approximate the targt function.
  llong PrefixSum(MultiplicitiveFunction mf,
                  const MultiplicitiveCombination& mc,
                  int mod);

 private:
  // Calculates \sum{ f(p) : 0 < p <= M and p is a prime } for some M.
  void GetSumOverPrimes(SumFunction mf_sum, int mod);

  void GetSumOverPrimes(const MultiplicitiveCombination& mc,
                        int mod);

  // Calculates the sum of all f(m) where:
  // - m is a composite number, and
  // - m is a multiple of X (not including given X), and
  // - X < m < N
  // Parameters:
  // - X: the given number to find multiples
  // - fx: the function value f(X)
  // - pIndex: The index of max prime factors of X (-1 if X = 1)
  // - mod: 0 if no mod is required, i.e. return real value
  llong GetSumOverMultiples(MultiplicitiveFunction mf,
                            llong X, llong fx, int pIndex,
                            int mod) const;

  llong getSumP(llong k) const {
    return k <= sqrt_n_ ? sum_p[k] : sum2_p[N/k];
  }

  void update(llong& a, llong index, int pIndex, int mod) {
    int p = primes[pIndex];
    a -= (getSumP(index/p) - getSumP(p-1)) * fp[pIndex];
    if (mod != 0) a %= mod;
  }

  llong N;
  int sqrt_n_;
  std::vector<int> primes;

  // Temporary varibles used during the calculation.

  // fp[i] is the f value for primes[i]. caching for better preformance
  std::unique_ptr<llong[]> fp;

  // sum_p[k] = \sum{ f(p) : 0 < p <= k && p is a prime }
  // sum2_p[k] = \sum{ f(p) : 0 < p <= N/k && p is a prime }
  // To simplify the calculation, we define 1 to be a prime too
  std::unique_ptr<llong[]> sum_p;
  std::unique_ptr<llong[]> sum2_p;
};

MultiplicitiveSum::MultiplicitiveSum(llong n) : N(n) {
  sqrt_n_ = int(sqrt(n));

  auto is_prime = std::make_unique<bool[]>(sqrt_n_+1);
  memset(is_prime.get(), true, sizeof(bool)*(sqrt_n_+1));

  for (int i = 2; i <= sqrt_n_; i++) if (is_prime[i]) {
    primes.push_back(i);
    for (int j = i+i; j <= sqrt_n_; j += i) is_prime[j] = false;
  }

  fp = std::make_unique<llong[]>(primes.size());
  sum_p = std::make_unique<llong[]>(sqrt_n_+1);
  sum2_p = std::make_unique<llong[]>(sqrt_n_+1);
}

void MultiplicitiveSum::GetSumOverPrimes(
    SumFunction mf_sum, int mod) {
  for (int i = 0; i < primes.size(); i++)
    fp[i] = mf_sum(primes[i], mod) - mf_sum(primes[i]-1, mod);
  for (int i = 0; i <= sqrt_n_; i++) sum_p[i] = mf_sum(i, mod);
  for (int i = 1; i <= sqrt_n_; i++) sum2_p[i] = mf_sum(N/i, mod);

  int pIndex = 0;
  for (llong p : primes) {
    // After each stage, sum_p and sum2_p contains all numbers that are:
    // - either a prime number, or
    // - a composite number whose smallest prime factor is >= primes[i]
    llong minN = p*(p-1);
    for (int j = 1; j <= sqrt_n_ && N/j > minN; j++)
      update(sum2_p[j], N/j, pIndex, mod);
    for (int j = sqrt_n_; j >= 0 && j > minN; j--)
      update(sum_p[j], j, pIndex, mod);
    pIndex++;
  }
}

void MultiplicitiveSum::GetSumOverPrimes(
    const MultiplicitiveCombination& mc, int mod) {
  std::unique_ptr<llong[]> sum_p_tmp_ =
      std::make_unique<llong[]>(sqrt_n_+1);
  std::unique_ptr<llong[]> sum2_p_tmp_ =
      std::make_unique<llong[]>(sqrt_n_+1);

  llong array_size = sizeof(llong)*(sqrt_n_+1);

  memset(sum_p_tmp_.get(), 0, array_size);
  memset(sum2_p_tmp_.get(), 0, array_size);

  for (auto entry : mc.GetFunctions()) {
    GetSumOverPrimes(std::get<1>(entry), mod);
    llong coef = std::get<0>(entry);
    for (int i = 0; i <= sqrt_n_; i++) {
      sum_p_tmp_[i] += coef*sum_p[i];
      sum2_p_tmp_[i] += coef*sum2_p[i];
      if (mod != 0) {
        sum_p_tmp_[i] %= mod;
        sum2_p_tmp_[i] %= mod;
      }
    }
  }

  memcpy(sum_p.get(), sum_p_tmp_.get(), array_size);
  memcpy(sum2_p.get(), sum2_p_tmp_.get(), array_size);
}

llong MultiplicitiveSum::GetSumOverMultiples(MultiplicitiveFunction mf,
    llong X, llong fx, int pIndex, int mod) const {
  llong res = 0;

  // Categorize the multiples Y into two buckets and calculate them separately:
  // 1. Y has only 1 largest prime factor
  // 2. Y has >= 2 largest prime factors
  for (int i = pIndex+1; i < primes.size(); i++) {
    int p = primes[i];
    if (1LL*p*p > X) break;

    llong nextX = X;
    for (int j = 1; nextX >= p; j++) {
      nextX /= p;
      llong nextF = fx * mf(p, j, mod);
      if (mod != 0) nextF %= mod;
      if (nextX > p) {
        // Summarize numbers in 1st category
        llong sumP = getSumP(nextX) - getSumP(p);
        res += nextF * sumP;
        if (mod != 0) res %= mod;
        res += GetSumOverMultiples(mf, nextX, nextF, i, mod);
      }
      // Summarize numbers in 2nd category
      if (j >= 2) res += nextF;
    }
  }
  return mod == 0 ? res : (res % mod + mod) % mod;
}

llong MultiplicitiveSum::PrefixSum(MultiplicitiveFunction mf,
                                   SumFunction mf_sum, int mod) {
  MultiplicitiveCombination::Builder builder;

  builder.AddFunction(1, mf_sum);
  return PrefixSum(mf, builder.Build(), mod);
}

llong MultiplicitiveSum::PrefixSum(MultiplicitiveFunction mf,
                                   const MultiplicitiveCombination& mc,
                                   int mod) {
  GetSumOverPrimes(mc, mod);

  llong res = GetSumOverMultiples(mf, N, 1, -1, mod) + getSumP(N) - getSumP(1) + 1;
  if (mod != 0) res = (res+mod) % mod;
  return res;
}

}  // namespace algo

#endif  // ALGO_MULTIPLICITIVE_SUM_H_
