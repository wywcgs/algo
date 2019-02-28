#ifndef ALGO_NUMBERS_H_
#define ALGO_NUMBERS_H_

#include <cassert>
#include <functional>
#include <memory>

#include "defs.h"
#include "modular.h"

namespace algo {

class Numbers {
 public:
  // Parameters:
  // - n: the max order to be processed
  // - modp: the mod, which must be a prime
  Numbers(int n_, int modp_) : Numbers(n_, modp_, false) {}
  Numbers(int n_, int mopp_, bool bernoulli_);

  int inv(int k) const { return inverse[k]; }
  // B_plus(k) = (-1)^k * B_minus(k)
  int B_plus(int k) const { return b_plus[k]; }
  int C(int n, int k) const {
    if (n < 0 || k < 0) return 0;
    else if (k == 0 || n == k) return 1;
    else return 1LL*factors[n]*factors_inv[k]%modp * factors_inv[n-k]%modp;
  }

  // Get \sum_{i=1}^n i^k mod P, where p is the prime in constructor
  int PrefixSum(llong n, int k) const;

 private:
  void GetInverse();
  void GetFactors();
  void GetFactorsInverse();
  void GetBernoulliPlus();

  std::unique_ptr<int[]> inverse;
  std::unique_ptr<int[]> factors;
  std::unique_ptr<int[]> factors_inv;
  std::unique_ptr<int[]> b_plus;

  int n;
  int modp;
  bool bernoulli;
};

Numbers::Numbers(int n_, int modp_, bool bernoulli_)
    : n(n_), modp(modp_), bernoulli(bernoulli_) {
  GetInverse();
  GetFactors();
  GetFactorsInverse();
  if (bernoulli) {
    GetBernoulliPlus();
  }
}

void Numbers::GetInverse() {
  inverse = std::make_unique<int[]>(n+2);
  inverse[0] = 0;
  inverse[1] = 1;
  for (int i = 2; i <= n+1; i++)
    inverse[i] = 1LL*(modp-modp/i)*inverse[modp%i]%modp;
}

void Numbers::GetFactors() {
  factors = std::make_unique<int[]>(n+2);
  factors[0] = 1;
  for (int i = 1; i <= n+1; i++)
    factors[i] = 1LL*i*factors[i-1]%modp;
}

void Numbers::GetFactorsInverse() {
  factors_inv = std::make_unique<int[]>(n+2);
  factors_inv[0] = 1;
  for (int i = 1; i <= n+1; i++)
    factors_inv[i] = 1LL*inverse[i]*factors_inv[i-1]%modp;
}

void Numbers::GetBernoulliPlus() {
  b_plus = std::make_unique<int[]>(n+1);
  b_plus[0] = 1;
  for (int i = 1; i <= n; i++) {
    llong b = 0;
    for (int j = 0; j < i; j++)
      b = (b+1LL*C(i+1, j)*b_plus[j])%modp;
    b_plus[i] = (modp+1-b*inverse[i+1]%modp)%modp;
  }
}

int Numbers::PrefixSum(llong n, int k) const {
  assert(bernoulli);

  n %= modp;

  llong r = 0;
  for (int i = 1; i <= k+1; i++) {
    llong c = 1LL*C(k+1, i)*B_plus(k+1-i)%modp;
    r = (r+c*powR(n, i, modp))%modp;
  }
  r = r*inv(k+1)%modp;
  return r;
}

}  // namespace algo

#endif  // ALGO_NUMBERS_H_
