#ifndef ALGO_PRIMES_H_
#define ALGO_PRIMES_H_

#include <algorithm>
#include <cassert>

#include "defs.h"
#include "modular.h"

namespace algo {
namespace {

// For given prime p, returns a number c which satisfies
// c^2 = -1 (mod p). Returns -1 if it doesn't exist.
int GetQuadraticNonresiduePrime(int p)
{
  if (p%4 != 1) return -1;

  for (int i = 2; i < p; i++) {
    int r = powR(i, (p-1)/4, p);
    if (1LL*r*r%p == p-1) return r;
  }
  return -1;
}

}  // namespace

// For given prime p, returns a pair of positive integers (a, b) which
// satisfies a^2 + b^2 = p and a <= b. And a = b holds only when p = 2.
PII TwoSquareSumDecompositionPrime(int p)
{
  if (p == 2) return PII(1, 1);
  if (p%4 != 1) return PII(-1, -1);

  int x = GetQuadraticNonresiduePrime(p);
  assert(1LL*x*x%p == p-1);

  if (1LL*x*x+1 == p) return PII(1, x);
  int a = p, b = x;
  while (true) {
    a %= b;
    std::swap(a, b);
    if (1LL*b*b > p) continue;
    a %= b;
    assert(1LL*b*b + 1LL*a*a == p);
    return PII(a, b);
  }
}

}  // namespace algo

#endif  // ALGO_PRIMES_H_
