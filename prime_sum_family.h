#ifndef ALGO_PRIME_SUM_FAMILY_H_
#define ALGO_PRIME_SUM_FAMILY_H_

#include <algorithm>
#include <cassert>
#include <cmath>
#include <vector>

#include "defs.h"

namespace algo {

template <typename T>
using PartialSumFunctionT = std::function<T(llong, int)>;

template <typename T>
class PrimeSumFamily {
 public:
  PrimeSumFamily(llong n, int m);

  // Calculates families of \sum{ f(p) : 0 < p <= N, p is a prime and p = a (mod M) }
  // for some N and M.
  // - mf_sum: the prefix sum function of an completely multiplicitive
  //     function f over number k*M+a. I.e., mf_sum(n, a) = (\sum_{i = a (mod M)} f(i)).
  void GetSumOverPrimes(PartialSumFunctionT<T> mf_sum);

  // Gets the \sum{ f(p) : 0 < p <= k, p is a prime and p = a (mod M) }.
  // only works for 0 <= a < M and gcd(a, M) = 1
  T GetSum(int a, llong k) const {
    int id = p2id[a];
    assert(id >= 0);
    return k <= SG_N ? sum[k][id] : sum2[N/k][id];
  }

 private:
  void update(std::vector<T>& sumv, llong index, int pIndex) {
    llong p = primes[pIndex];

    std::vector<T> psum(id2p.size(), T(0));
    for (int i = 0; i < id2p.size(); i++)
      psum[i] = GetSum(id2p[i], index/p) - GetSum(id2p[i], p-1);

    for (int a : id2p) {
      int prev = p2id[a], cur = p2id[a*p%M];
      sumv[prev] -= psum[cur] * fp[pIndex];
    }
  }

  llong N;
  // SG_N = floor(N^{1/2})
  int SG_N;
  // Sums over Z/Zm field
  int M;

  std::vector<int> primes;

  // fp[i] is the f value for primes[i]. caching for better preformance
  std::vector<T> fp;

  // sum[i][k] is the prefix sums of k-th number in id2p
  std::vector<std::vector<T>> sum;
  std::vector<std::vector<T>> sum2;

  // mappings
  // id of the co-prime number to co-prime number.
  std::vector<int> id2p;
  // co-prime number to its id. Only defined in number coprime to m.
  std::vector<int> p2id;
};

template <typename T>
PrimeSumFamily<T>::PrimeSumFamily(llong n, int m) : N(n), M(m) {
  SG_N = int(sqrt(n));

  auto is_prime = std::vector<bool>(SG_N+1, true);
  for (int i = 2; i <= SG_N; i++) if (is_prime[i]) {
    primes.push_back(i);
    for (llong j = 1LL*i*i; j <= SG_N; j += i) is_prime[j] = false;
  }

  p2id = std::vector<int>(M, -1);
  for (int i = 0; i < M; i++) {
    if (std::__gcd(i, M) != 1) continue;
    int id = id2p.size();
    id2p.push_back(i);
    p2id[i] = id;
  }

  fp = std::vector<T>(primes.size());
  sum = std::vector<std::vector<T>>(SG_N+1, std::vector<T>(id2p.size(), T(0)));
  sum2 = std::vector<std::vector<T>>(SG_N+1, std::vector<T>(id2p.size(), T(0)));
}

template <typename T>
void PrimeSumFamily<T>::GetSumOverPrimes(PartialSumFunctionT<T> mf_sum) {
  for (int i = 0; i < primes.size(); i++) {
    int p = primes[i];
    if (std::__gcd(p, M) != 1) continue;
    fp[i] = mf_sum(p, p%M) - mf_sum(p-1, p%M);
  }
  for (int a : id2p) {
    int id = p2id[a];
    for (int i = 0; i <= SG_N; i++) sum[i][id] = mf_sum(i, a);
    for (int i = 1; i <= SG_N; i++) sum2[i][id] = mf_sum(N/i, a);
  }

  for (int pIndex = 0; pIndex < primes.size(); pIndex++) {
    llong p = primes[pIndex];
    if (std::__gcd<int>(p, M) != 1) continue;

    llong minN = p*(p-1);
    for (int j = 1; j <= SG_N && N/j > minN; j++)
      update(sum2[j], N/j, pIndex);
    for (int j = SG_N; j >= 0 && j > minN; j--)
      update(sum[j], j, pIndex);
  }
}

}  // namespace algo

#endif  // ALGO_PRIME_SUM_FAMILY_H_
