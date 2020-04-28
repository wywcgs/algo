#ifndef ALGO_PRIMES_H_
#define ALGO_PRIMES_H_

#include <algorithm>
#include <vector>
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

typedef __int128 i128;

llong powK(llong a, llong b, llong p)
{
  if (b == 0) return 1%p;
  if (b == 1) return a%p;
  i128 x = powK(a, b/2, p);
  x = x*x%p;
  if (b%2 == 1) x = x*a%p;
  return x;
}

// Run Rabin-Miller testing once against given prime "a"
bool RmOnce(llong n, int a)
{
  llong d = n - 1;
  int s = 1;
  while (d%2 == 0) {
    d /= 2;
    s++;
  }

  llong x = powK(a, d, n);
  if (x == 1 || x == n-1) return true;

  for (int i = 0; i < s-1; i++) {
    x = (i128)x*x%n;
    if (x == n-1) return true;
  }
  return false;
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

// Rabin Miller prime testing
bool RabinMiller(llong n)
{
  const std::vector<int> a = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
  if (n <= a.back()) return std::find(a.begin(), a.end(), n) != a.end();

  for (int i : a) if (!RmOnce(n, i)) return false;
  return true;
}

}  // namespace algo
