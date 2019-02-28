#ifndef ALGO_MULTIPLICITIVE_SUM2_H_
#define ALGO_MULTIPLICITIVE_SUM2_H_

#include <algorithm>
#include <cmath>
#include <memory>
#include <vector>

#include "defs.h"
#include "division_enumerator.h"

namespace algo {

template<typename T>
class MultiplicitiveSum2 {
 public:
  MultiplicitiveSum2(llong n) : N(n) {
    BF_N = llong(pow(N, 2.0/3));
    SG_N = int(pow(N, 1.0/2));

    sum1 = std::make_unique<T[]>(SG_N+1);
    sum2 = std::make_unique<T[]>(SG_N+1);
    smallSum = std::make_unique<T[]>(BF_N+1);
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
  void Calculate(MultiplicitiveFunctionT<T> f,
      NtFunctionT<T> gPrefixSum, NtFunctionT<T> rPrefixSum);

  // Gets the prefix sum of index m. Must be called after method Calculate.
  T GetPrefixSum(llong m) const { return m <= BF_N ? smallSum[m] : sum2[N/m]; }
 
 private:
  void CalculateSmallSums(MultiplicitiveFunctionT<T> f);
 
  T CalculateIndex(NtFunctionT<T> gPrefixSum,
      NtFunctionT<T> rPrefixSum, llong m) const;
 
  llong N;
  // BF_N = floor(N^{2/3})
  llong BF_N;
  // SG_N = floor(N^{1/2})
  int SG_N;

  std::unique_ptr<T[]> sum1;
  std::unique_ptr<T[]> sum2;
  std::unique_ptr<T[]> smallSum;
};

template <typename T>
void MultiplicitiveSum2<T>::CalculateSmallSums(MultiplicitiveFunctionT<T> f) {
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
  smallSum[0] = T(0);
  smallSum[1] = T(1);
  for (llong i = 2; i <= BF_N; i++) {
    if (prime[i]) {
      smallSum[i] = f(i, 1);
      continue;
    }

    int pc = 0, p = pFactor[i];
    llong m = i;
    for (; m % p == 0; m /= p) pc++;
    smallSum[i] = smallSum[m] * f(p, pc);
  }
  
  // Get the prefix sums.
  for (int i = 1; i <= BF_N; i++)
    smallSum[i] += smallSum[i-1];
}

template <typename T>
void MultiplicitiveSum2<T>::Calculate(MultiplicitiveFunctionT<T> f,
    NtFunctionT<T> gPrefixSum, NtFunctionT<T> rPrefixSum) {
  CalculateSmallSums(f);

  for (int i = 1; i <= SG_N; i++)
    sum1[i] = CalculateIndex(gPrefixSum, rPrefixSum, i);

  for (int i = SG_N; i >= 1; i--)
    sum2[i] = CalculateIndex(gPrefixSum, rPrefixSum, N/i);
}

template <typename T>
T MultiplicitiveSum2<T>::CalculateIndex(
    NtFunctionT<T> gPrefixSum, NtFunctionT<T> rPrefixSum, llong m) const {
  if (m <= BF_N) return smallSum[m];

  DivisionEnumerator de;
  T res = rPrefixSum(m);

  de.Do(m, [this, &res, &rPrefixSum, &gPrefixSum]
      (llong coef, llong begin, llong end) {
    T coefSum = gPrefixSum(end) - gPrefixSum(begin);
    T fSum = GetPrefixSum(coef);
    res -= coefSum * fSum;
  });

  return res;
}

}  // namespace algo

#endif  // ALGO_MULTIPLICITIVE_SUM2_H_
