#ifndef ALGO_INTEGRAL_H_
#define ALGO_INTEGRAL_H_

#include <functional>
#include <cassert>

namespace algo {

template<typename T>
T simpson(std::function<T(T)> f, T lo, T hi, int n)
{
  assert(n > 0);
  if (n%2 != 0) n++;

  T dx = (hi-lo)/n;
  T res = f(lo) + f(hi);
  for (int i = 2; i < n; i += 2) res += 2*f(lo+i*dx);
  for (int i = 1; i < n; i += 2) res += 4*f(lo+i*dx);
  res *= dx/3;

  return res;
}

}  // namespace algo

#endif  // ALGO_INTEGRAL_H_
