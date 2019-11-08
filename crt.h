#ifndef ALGO_CRT_H_
#define ALGO_CRT_H_

#include <algorithm>

#include "gcd.h"

namespace algo {

template <typename T>
class CRT {
 public:
  CRT(T a, T n): a_(a), n_(n) {}
  CRT() : CRT(0, 1) {}

  bool Merge(T a, T n);
  T GetSolution() const { return a_; }
  T GetModule() const { return n_; }

 private:
  T a_;
  T n_;
};

template <typename T>
bool CRT<T>::Merge(T a2, T n2) {
  T g = std::__gcd(n_, n2);
  T c = a2-a_;
  if (c % g != 0) return false;

  T x1, x2;
  ExtendGcd(n_/g, n2/g, c/g, x1, x2);
  a_ += n_*x1;
  n_ *= n2/g;
  return true;
}

}  // namespace algo

#endif  // ALGO_CRT_H_
