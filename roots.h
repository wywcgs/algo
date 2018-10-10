#ifndef ALGO_ROOTS_H_
#define ALGO_ROOTS_H_

#include <vector>

#include "modular.h"

namespace algo {

namespace {

int GetPower(int p, int e)
{
  int res = 1;
  for (int i = 0; i < e; i++) res *= p;
  return res;

};

}  // namespace

// Gets the smallest positive number i which satisfies r^i = 1 (mod p^e)
int GetOrder(int r, int p, int e, const std::vector<int>& factors_of_phi_p)
{
  int pmod = GetPower(p, e), phi = pmod/p*(p-1);
  
  int order = phi;
  while(order%p == 0 && powR(r, order/p, pmod) == 1) order /= p;
  for (int factor : factors_of_phi_p)
    while (order%factor == 0 && powR(r, order/factor, pmod) == 1)
      order /= factor;
  return order;
}

// Gets the primitive root of p^e. The prime factors of of (p-1) must be given.
// Restrictions:
// - p must be an odd prime number
// - p^e must be smaller than 32 bit integer max
int PrimitiveRoot(int p, int e, const std::vector<int>& factors_of_phi_p)
{
  int phi = GetPower(p, e)/p*(p-1);

  for (int i = 2; true; i++)
    if (GetOrder(i, p, e, factors_of_phi_p) == phi)
      return i;
}

int PrimitiveRoot(int p, const std::vector<int>& factors_of_phi_p)
{
  return PrimitiveRoot(p, 1, factors_of_phi_p);
}

}  // namespace algo

#endif  // ALGO_ROOTS_H_
