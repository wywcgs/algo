#ifndef ALGO_ROOTS_H_
#define ALGO_ROOTS_H_

#include <functional>
#include <vector>

#include "modular.h"

namespace algo {

namespace {

llong GetPower(llong p, int e)
{
  llong res = 1;
  for (int i = 0; i < e; i++) res *= p;
  return res;

};

}  // namespace

// Gets the smallest positive number i which satisfies r^i = 1 (mod p^e)
llong GetOrder(llong r, llong p, int e, const std::vector<int>& factors_of_phi_p)
{
  llong pmod = GetPower(p, e), phi = pmod/p*(p-1);
  
  llong order = phi;
  while(order%p == 0 && powR64(r, order/p, pmod) == 1) order /= p;
  for (int factor : factors_of_phi_p)
    while (order%factor == 0 && powR64(r, order/factor, pmod) == 1)
      order /= factor;
  return order;
}

// Gets the primitive root of p^e. The prime factors of of (p-1) must be given.
// Restrictions:
// - p must be an odd prime number
// - p^e must be smaller than 10^18
int PrimitiveRoot(llong p, int e, const std::vector<int>& factors_of_phi_p)
{
  llong phi = GetPower(p, e)/p*(p-1);

  for (int i = 2; true; i++)
    if (GetOrder(i, p, e, factors_of_phi_p) == phi)
      return i;
}

int PrimitiveRoot(llong p, const std::vector<int>& factors_of_phi_p)
{
  return PrimitiveRoot(p, 1, factors_of_phi_p);
}

// Find one solution to x^2 = n (mod p) using Tonelli-Shanks algorithm
// Restrictions:
// - the root must exist, i.e., n^{(p-1)/2} = 1 (mod p)
// - p must be a prime number
// - p must be smaller than 10^18
// If no solutions, the method will return -1.
llong SquareRoot(llong n, llong p, const std::vector<int>& factors_of_phi_p)
{
  if (n == 0) return 0;

  // Find p-1 = Q*2^S
  llong s = 0, q = p-1;
  while (q%2 == 0) {
    q /= 2;
    s++;
  }

  llong z = PrimitiveRoot(p, factors_of_phi_p);
  llong m = s, c = powR64(z, q, p), t = powR64(n, q, p), r = powR64(n, (q+1)/2, p);
  
  while (t != 1) {
    //printf("n = %lld, p = %lld, s = %lld, q = %lld, t = %lld, c = %lld, r = %lld, m = %lld\n", n, p, s, q, t, c, r, m);
    int i = 1;
    llong ct = multiply64(t, t, p);
    for (; i < m && ct != 1; i++) ct = multiply64(ct, ct, p);
    if (i == m) return -1;

    llong b = powR64(c, 1LL<<(m-i-1), p);
    m = i;
    c = multiply64(b, b, p);
    t = multiply64(t, multiply64(b, b, p), p);
    r = multiply64(r, b, p);
  }

  return r;
}

// Given a solution r to f(r) = 0 (mod p^e), find the solution r' to
// f(r') = 0 (mod p^{e+1}) using Hensel's lemma.
// 
// Parameters:
// - f(r, p) calculates f(r) mod p
// - derive_f(r, p) calculates f'(r) mod p 
//
// Restrictions:
// - r must be a solution to f(r) = 0 (mod p^e), i.e., f(r, p^e) = 0
// - p must be a prime number
// - p^{e+1} must be smaller than 10^18
llong LiftSolution(llong r, llong p, int e,
                   std::function<llong(llong, llong)> f,
                   std::function<llong(llong, llong)> derive_f)
{
  llong pnow = GetPower(p, e), pmod = pnow*p;
  llong fr = f(r, pmod), dfr = derive_f(r, pmod) % p;

  if (dfr == 0) return -1;
  llong t = multiply64(p-fr/pnow, inverse(dfr, p), p);
  
  return r + t*pnow;
}

}  // namespace algo

#endif  // ALGO_ROOTS_H_
