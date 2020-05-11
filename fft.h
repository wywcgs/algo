#ifndef ALGO_FFT_H_
#define ALGO_FFT_H_

#include <iostream>
#include <complex>
#include <cmath>
#include <vector>

#include "defs.h"
#include "mod_num.h"

namespace algo {
namespace {

const double PI = acos(-1);

template <typename T>
class ArraySlice {
 public:
  ArraySlice(std::vector<T>* pA, int psz, int pa = 1, int pb = 0)
      : A(pA), sz(psz), a(pa), b(pb) {
    assert(A != nullptr);
  }
  ArraySlice(std::vector<T>& pA, int psz) : ArraySlice(&pA, psz) {}
  ArraySlice(std::vector<T>& pA) : ArraySlice(&pA, pA.size()) {}

  // Subarray which contains even indexes of current one.
  ArraySlice Even() { return ArraySlice(A, sz/2, 2*a, b); }
  // Subarray which contains odd indexes of current one.
  ArraySlice Odd() { return ArraySlice(A, sz-sz/2, 2*a, b+a); }
  // Subarray of the first half. Must apply on even length array slices.
  ArraySlice FirstHalf() { return ArraySlice(A, sz/2, a, b); }
  // Subarray of the second half. Must apply on even length array slices.
  ArraySlice SecondHalf() { return ArraySlice(A, sz/2, a, b+sz/2); }

  int size() const { return sz; }
  T& operator[](int n) { assert(a*n+b < A->size()); return (*A)[a*n+b]; }
  const T operator[](int n) const { assert(a*n+b < A->size()); return (*A)[a*n+b]; }

 private:
  std::vector<T>* A; // not owned
  // i-th element is A[a*i+b]
  int sz;
  int a;
  int b;
};

template <typename T>
void FftInternal(ArraySlice<T> a, ArraySlice<T> c, ArraySlice<T> w)
{
  int n = w.size();
  if (n == 1) {
    c[0] = a[0];
    return;
  }

  int hn = n/2;
  FftInternal(a.Even(), c.FirstHalf(), w.Even());
  FftInternal(a.Odd(), c.SecondHalf(), w.Even());
  for (int i = 0; i < hn; i++) {
    T a1 = c[i], b1 = c[i+hn];
    c[i] = a1 + w[i] * b1;
    c[i+hn] = a1 - w[i] * b1;
  }
}

}  // namespace


std::vector<llong> fft(const std::vector<int>& pa, const std::vector<int>& pb)
{
  typedef std::complex<double> ftype;

  int len = pa.size() + pb.size() - 1, n = 1;
  for (; n < len; n *= 2)
    ;

  std::vector<ftype> a(n, ftype()), b(n, ftype());
  for (int i = 0; i < pa.size(); i++) a[i] = ftype(pa[i]);
  for (int i = 0; i < pb.size(); i++) b[i] = ftype(pb[i]);

  auto root = std::vector<ftype>(n);
  for (int i = 0; i < n; i++) root[i] = std::polar(1.0, 2*PI/n*i);
  auto r = std::vector<ftype>(n);
  auto rt = std::vector<ftype>(n);

  double sqrtn = sqrt(n);
  //FftInternal(as, rts, ws);
  FftInternal<ftype>(a, rt, root);
  for (int i = 0; i < n; i++) r[i] = rt[i] / sqrtn;
  FftInternal<ftype>(b, rt, root);
  for (int i = 0; i < n; i++) r[i] *= rt[i] / sqrtn;

  // Get inverse of w
  reverse(root.begin()+1, root.end());
  FftInternal<ftype>(r, rt, root);

  std::vector<llong> res(n);
  for (int i = 0; i < n; i++) {
    res[i] = llong(round(rt[i].real()));
//    printf("%d %lld\n", i, res[i]);
  }
  return res;
}

}  // namespae algo

#endif  // ALGO_FFT_H_
