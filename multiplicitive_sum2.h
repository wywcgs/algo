#ifndef ALGO_MULTIPLICITIVE_SUM2_H_
#define ALGO_MULTIPLICITIVE_SUM2_H_

#include <algorithm>
#include <cmath>
#include <memory>
#include <vector>

#include "defs.h"
#include "division_enumerator.h"

namespace algo {

class MultiplicitiveSum2 {
 public:
  MultiplicitiveSum2(llong n) : N(n) {
    BF_N = llong(pow(N, 2.0/3));
    SG_N = int(pow(N, 1.0/2));
    
    sum1 = std::make_unique<llong[]>(SG_N+1);
    sum2 = std::make_unique<llong[]>(SG_N+1);
    smallSum = std::make_unique<llong[]>(BF_N+1);
  }

  // Calculates the prefix sum of given multiplicitive function in some fixed prefixes.
  // 
  // Input:
  // - f: the original function which you want to get prefix sum. You only need to
  //      define it on all number prime numbers and their powers.
  // - gPrefixSum: the prefix sum function of function g, which you apply Dirichlet
  //               convolution with original function.
  // - rPrefixSum: the prefix sum function for Dirichlet convolution r. So the following
  //               equation must be held for all params: f * g = r. Here '*' is the
  //               Dirichlet convolution.
  //
  // Output: 
  // - It can calculate $2*sqrt(N)$ prefix sums alltogether. A prefix index m can be
  //   calculated if there is another number m' which satifies floor(N/m') = m.
  //   The result can be accessed via "GetPrefixSum" method. 
  void Calculate(MultiplicitiveFunction f,
      NtFunction gPrefixSum, NtFunction rPrefixSum, int modp);

  // Gets the prefix sum of index m. Must be called after method Calculate.
  llong GetPrefixSum(llong m) const { return m <= BF_N ? smallSum[m] : sum2[N/m]; }
 
 private:
  void CalculateSmallSums(MultiplicitiveFunction f, int modp);
 
  llong CalculateIndex(NtFunction gPrefixSum,
      NtFunction rPrefixSum, llong m, int modp) const;
 
  llong N;
  // BF_N = floor(N^{2/3})
  llong BF_N;
  // SG_N = floor(N^{1/2})
  int SG_N;

  std::unique_ptr<llong[]> sum1;
  std::unique_ptr<llong[]> sum2;
  std::unique_ptr<llong[]> smallSum;
};

void MultiplicitiveSum2::CalculateSmallSums(MultiplicitiveFunction f, int modp) {
  std::vector<bool> prime(BF_N+1, true);
  std::vector<int> pFactor(BF_N+1, 0);

  for (llong i = 2; i <= BF_N; i++) if (prime[i]) {
    llong start = i*i;
    if (start > BF_N) break;
    for (llong j = start; j <= BF_N; j += i) {
      prime[j] = false;
      pFactor[j] = i;
    }
  }

  // Calculate the value of function f in smallSum.
  smallSum[0] = 0;
  smallSum[1] = 1;
  for (llong i = 2; i <= BF_N; i++) {
    if (prime[i]) {
      smallSum[i] = f(i, 1, modp);
      continue;
    }

    int pc = 0, p = pFactor[i];
    llong m = i;
    for (; m % p == 0; m /= p) pc++;
    smallSum[i] = smallSum[m] * f(p, pc, modp);
    if (modp != 0) smallSum[i] %= modp;
  }
  
  // Get the prefix sums.
  for (int i = 1; i <= BF_N; i++) {
    smallSum[i] += smallSum[i-1];
    if (modp != 0) smallSum[i] %= modp;
  }
}

void MultiplicitiveSum2::Calculate(MultiplicitiveFunction f,
    NtFunction gPrefixSum, NtFunction rPrefixSum, int modp) {
  CalculateSmallSums(f, modp);

  for (int i = 1; i <= SG_N; i++)
    sum1[i] = CalculateIndex(gPrefixSum, rPrefixSum, i, modp);

  for (int i = SG_N; i >= 1; i--)
    sum2[i] = CalculateIndex(gPrefixSum, rPrefixSum, N/i, modp);
}

llong MultiplicitiveSum2::CalculateIndex(
    NtFunction gPrefixSum, NtFunction rPrefixSum, llong m, int modp) const {
  if (m <= BF_N) return smallSum[m];

  DivisionEnumerator de;
  llong res = rPrefixSum(m, modp);

  de.Do(m, [this, &res, &rPrefixSum, &gPrefixSum, &modp]
      (llong coef, llong begin, llong end) {
    llong coefSum = gPrefixSum(end, modp) - gPrefixSum(begin, modp);
    llong fSum = GetPrefixSum(coef);
    res -= coefSum * fSum;
    if (modp != 0) res %= modp;
  });

  return res;
}

}  // namespace algo

#endif  // ALGO_MULTIPLICITIVE_SUM2_H_
