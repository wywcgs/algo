#ifndef ALGO_MODULAR_H_
#define ALGO_MODULAR_H_

#include <algorithm>
#include <cassert>

#include "defs.h"
#include "gcd.h"

namespace algo {

// Multiply big numbers.
// Restrictions:
// - 0 <= c < 2^60 (~1e18)
llong multiply64(llong a, llong b, llong c)
{
  if (a == 0 || b == 0) return 0;
  if (a == 1 || b == 1) return a*b%c;

  a %= c;
  b %= c;

  const llong MAX_V = 1LL<<61;
  llong res = 0, d = MAX_V/c;
  while (b > 0) {
    res = (res+b%d*a) % c;
    a = (a*d)%c;
    b /= d;
  }

  return res;
}

llong powR64(llong a, llong n, llong p)
{
  if (n == 0) return 1;
  if (n == 1) return a%p;

  llong r = powR64(a, n>>1, p);
  r = multiply64(r, r, p);
  if (n&1) r = multiply64(r, a, p);
  return r;
}

template <typename T>
int powR(int a, T n, int p)
{
  if (n == 0) return 1;
  if (n == 1) return a%p;

  llong r = powR(a, n>>1, p);
  r = r*r%p;
  if (n&1) r = r*a%p;
  return r;
}
  
template <typename V, typename T>
V powR(V a, T n)
{
  if (n == 0) return 1;
  if (n == 1) return a;

  V r = powR(a, n>>1);
  r = r*r;
  if (n&1) r = r*a;
  return r;
}

// Gets number b in (0, p) where a*b = 1 (mod p).
int inverse(int a, int p)
{
  // a*x + p*y = 1
  llong x, y;
  assert(ExtendGcd<llong>(a, p, 1, x, y));
  return x;
}

// Gets number b in (0, p) where a*b = 1 (mod p).
// p must be a prime to work
llong inverse(llong a, llong p)
{
  return powR64(a, p-2, p);
}

}  // namespace algo

#endif  // ALGO_MODULAR_H_
