#ifndef ALGO_MOD_NUM_H
#define ALGO_MOD_NUM_H

#include <functional>

#include "defs.h"
#include "modular.h"

namespace algo {

template <int P>
class ModNum {
 public:
  ModNum() : n(0) {}
  ModNum(llong n_) {
    if (0 <= n_ && n_ < P) n = n_;
    else if (n_ >= P && n_ < P+P) n = n_-P;
    else if (n_ < 0 && n_+P >= 0) n = n_+P;
    else {
      int m = n_%P;
      n = (m >= 0) ? m : P-m;
    }
  }

  operator int() const { return n; }

  // Arithmic operations
  // TODO: support division operators
  friend ModNum<P> operator +(const ModNum<P>& lhs, const ModNum<P>& rhs) { return ModNum<P>(lhs.n+rhs.n); }
  friend ModNum<P> operator +(int lhs, const ModNum<P>& rhs) { return ModNum<P>(lhs+rhs.n); }
  friend ModNum<P> operator +(const ModNum<P>& lhs, int rhs) { return ModNum<P>(lhs.n+rhs); }
  friend ModNum<P> operator +(llong lhs, const ModNum<P>& rhs) { return ModNum<P>(lhs+rhs.n); }
  friend ModNum<P> operator +(const ModNum<P>& lhs, llong rhs) { return ModNum<P>(lhs.n+rhs); }
  
  friend ModNum<P> operator -(const ModNum<P>& lhs, const ModNum<P>& rhs) { return ModNum<P>(lhs.n-rhs.n); }
  friend ModNum<P> operator -(int lhs, const ModNum<P>& rhs) { return ModNum<P>(lhs-rhs.n); }
  friend ModNum<P> operator -(const ModNum<P>& lhs, int rhs) { return ModNum<P>(lhs.n-rhs); }
  friend ModNum<P> operator -(llong lhs, const ModNum<P>& rhs) { return ModNum<P>(lhs-rhs.n); }
  friend ModNum<P> operator -(const ModNum<P>& lhs, llong rhs) { return ModNum<P>(lhs.n-rhs); }
  
  friend ModNum<P> operator *(const ModNum<P>& lhs, const ModNum<P>& rhs) { return ModNum<P>(1LL*lhs.n*rhs.n); }
  friend ModNum<P> operator *(int lhs, const ModNum<P>& rhs) { return ModNum<P>(1LL*lhs*rhs.n); }
  friend ModNum<P> operator *(const ModNum<P>& lhs, int rhs) { return ModNum<P>(1LL*lhs.n*rhs); }
  friend ModNum<P> operator *(llong lhs, const ModNum<P>& rhs) { return ModNum<P>(lhs)*rhs; }
  friend ModNum<P> operator *(const ModNum<P>& lhs, llong rhs) { return lhs*ModNum<P>(rhs); }
  
  ModNum<P>& operator += (const ModNum<P>& b) { return *this = *this + b; }
  ModNum<P>& operator += (int b) { return *this = *this + b; }
  ModNum<P>& operator += (llong b) { return *this = *this + b; }

  ModNum<P>& operator -= (const ModNum<P>& b) { return *this = *this - b; }
  ModNum<P>& operator -= (int b) { return *this = *this - b; }
  ModNum<P>& operator -= (llong b) { return *this = *this - b; }
  
  ModNum<P>& operator *= (const ModNum<P>& b) { return *this = *this * b; }
  ModNum<P>& operator *= (int b) { return *this = *this * b; }
  ModNum<P>& operator *= (llong b) { return *this = *this * b; }
  
  ModNum<P> pow(llong m) const { return ModNum<P>(powR(n, m, P)); }

 private:
  int n;
};

}  // namespace algo

#endif  // ALGO_MOD_NUM_H
