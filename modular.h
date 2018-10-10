#ifndef ALGO_MODULAR_H_
#define ALGO_MODULAR_H_

namespace algo {

#include "defs.h"

int powR(int a, llong n, int p)
{
  if (n == 0) return 1;
  if (n == 1) return a%p;

  llong r = powR(a, n>>1, p);
  r = r*r%p;
  if (n&1) r = r*a%p;
  return r;
}

// Gets number b in (0, p) where a*b = 1 (mod p).
// p must be a prime to work
int inverse(int a, int p)
{
  return powR(a, p-2, p);
}

};

#endif  // ALGO_MODULAR_H_
