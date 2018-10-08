#ifndef ALGO_H_
#define ALGO_H_

namespace algo {

#include <algorithm>
 
template <typename T>
class CRT {
 public:
  CRT(T a, T n): a_(a), n_(n) {}
  CRT() : CRT(0, 1) {}
  bool Merge(T a, T n);
  T GetSolution() const { return a_; }
  T GetModule() const { return n_; }

 private:
  void ExtGcd(T n1, T n2, T c, T& x1, T& x2) const;

  T a_;
  T n_;
};

template <typename T>
void CRT<T>::ExtGcd(T n1, T n2, T c, T& x1, T& x2) const {
  if (n2 == 0) {
    x1 = c;
    x2 = 0;
  } else {
    T y1, y2;
    // n1*x1 + n2*x2 = c
    // ==> n2*x2 + (n1%n2 + k*n2)*x1 = c
    // ==> n2*(x2+k*x1) + n1%n2*x1 = c
    // ==> x1 = y2, x2 = y1 - k*y2
    ExtGcd(n2, n1%n2, c, y1, y2);
    x1 = y2;
    x2 = y1 - n1/n2*y2;

    T delta = x1/n2;
    x1 -= delta*n2;
    x2 += delta*n1;
    if (x1 < 0) {
      x1 += n2;
      x2 -= n1;
    }
  }
}

template <typename T>
bool CRT<T>::Merge(T a2, T n2) {
  T g = std::__gcd(n_, n2);
  T c = a2-a_;
  if (c % g != 0) return false;

  T x1, x2;
  ExtGcd(n_/g, n2/g, c/g, x1, x2);
  a_ += n_*x1;
  n_ *= n2/g;
  return true;
}

};

#endif  // ALGO_H_
