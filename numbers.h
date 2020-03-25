#ifndef ALGO_NUMBERS_H_
#define ALGO_NUMBERS_H_

#include <cassert>
#include <functional>
#include <memory>

#include "defs.h"
#include "modular.h"

namespace algo {

template <typename T>
class Numbers {
 public:
  // Parameters:
  // - n: the max order to be processed
  Numbers(int n_) : Numbers(n_, false) {}
  Numbers(int n_, bool bernoulli_);

  T inv(int k) const { return inversion[k]; }
  // B_plus(k) = (-1)^k * B_minus(k)
  T B_plus(int k) const { return b_plus[k]; }
  T C(int n, int k) const {
    if (n < 0 || k < 0 || n < k) return T(0);
    if (k == 0 || n == k) return T(1);
    return factors[n] * factors_inv[k] * factors_inv[n-k];
  }

  // Get \sum_{i=1}^n i^k mod P, where p is the prime in constructor
  T PrefixSum(llong n, int k) const;
  // Get the coefficient of the term "n^k" in the expression \sum_{i=1}^n i^m
  T PrefixSumCoef(int m, int k) const;

 private:
  void GetInverse();
  void GetFactors();
  void GetFactorsInverse();
  void GetBernoulliPlus();

  std::unique_ptr<T[]> inversion;
  std::unique_ptr<T[]> factors;
  std::unique_ptr<T[]> factors_inv;
  std::unique_ptr<T[]> b_plus;

  int n;
  bool bernoulli;
};

template<typename T>
Numbers<T>::Numbers(int n_, bool bernoulli_)
    : n(n_), bernoulli(bernoulli_) {
  GetInverse();
  GetFactors();
  GetFactorsInverse();
  if (bernoulli) {
    GetBernoulliPlus();
  }
}

template<typename T>
void Numbers<T>::GetInverse() {
  inversion = std::make_unique<T[]>(n+2);
  inversion[0] = T(0);
  inversion[1] = T(1);
  for (int i = 2; i <= n+1; i++)
    //inversion[i] = 1LL*(modp-modp/i)*inverse[modp%i]%modp;
    inversion[i] = 1/T(i);
}

template<typename T>
void Numbers<T>::GetFactors() {
  factors = std::make_unique<T[]>(n+2);
  factors[0] = T(1);
  for (int i = 1; i <= n+1; i++)
    factors[i] = i*factors[i-1];
}

template<typename T>
void Numbers<T>::GetFactorsInverse() {
  factors_inv = std::make_unique<T[]>(n+2);
  factors_inv[0] = T(1);
  for (int i = 1; i <= n+1; i++)
    factors_inv[i] = inversion[i]*factors_inv[i-1];
}

template<typename T>
void Numbers<T>::GetBernoulliPlus() {
  b_plus = std::make_unique<T[]>(n+1);
  b_plus[0] = 1;
  for (int i = 1; i <= n; i++) {
    T b = 0;
    for (int j = 0; j < i; j++)
      b += C(i+1, j)*b_plus[j];
    b_plus[i] = 1 - b*inversion[i+1];
  }
}

template<typename T>
T Numbers<T>::PrefixSum(llong n, int k) const {
  assert(bernoulli);

  T r = 0, pn = 1;
  for (int i = 1; i <= k+1; i++) {
    pn *= n;
    r += PrefixSumCoef(k, i)*pn;
  }
  return r;
}

template<typename T>
T Numbers<T>::PrefixSumCoef(int m, int k) const {
  assert(bernoulli);

  if (k <= 0) return T(0);
  return C(m+1, k)*B_plus(m+1-k)*inversion[m+1];
}

}  // namespace algo

#endif  // ALGO_NUMBERS_H_
