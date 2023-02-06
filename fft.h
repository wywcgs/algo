#ifndef ALGO_FFT_H_
#define ALGO_FFT_H_

#include <complex>
#include <cmath>
#include <tuple>
#include <vector>

#include "defs.h"
#include "modular.h"
#include "mod_num.h"

namespace algo {

enum FwhtOperator { XOR, AND, OR };

namespace {

typedef std::complex<ldouble> ftype;

const ldouble PI = acosl(-1);

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

template <typename T>
void FwhtInternal(std::vector<T>& a, FwhtOperator op, bool inverse)
{
  int n = a.size();
  for (int l = 1; l < n; l *= 2) {
    for (int i = 0; i < n; i += l+l) {
      for (int j = i; j < i+l; j++) {
        T &x = a[j], &y = a[j+l];
        switch (op) {
          case XOR:
            std::tie(x, y) = std::make_pair(x+y, x-y);
            break;
          case AND:
            std::tie(x, y) = inverse ? std::make_pair(x-y, y) : std::make_pair(x, x+y);
            break;
          case OR:
            std::tie(x, y) = inverse ? std::make_pair(x, x-y) : std::make_pair(x+y, y);
          default:
            // should never happen.
            assert(false);
        }
      }
    }
  }
  if (inverse && op == XOR) {
    for (T& i : a) i /= n;
  }
}

int FindFftSize(const std::vector<int>& pa, const std::vector<int>& pb)
{
  int len = pa.size() + pb.size() - 1, n = 1;
  for (; n < len; n *= 2)
    ;
  return n;
}

}  // namespace


std::vector<llong> fft(const std::vector<int>& pa, const std::vector<int>& pb)
{
  int n = FindFftSize(pa, pb);
  ldouble sqrtn = sqrt(n);
  std::vector<ftype> a(n, ftype()), b(n, ftype());
  for (int i = 0; i < pa.size(); i++) a[i] = ftype(pa[i]) / sqrtn;
  for (int i = 0; i < pb.size(); i++) b[i] = ftype(pb[i]) / sqrtn;

  auto root = std::vector<ftype>(n);
  for (int i = 0; i < n; i++) root[i] = std::polar(1.0L, 2*PI/n*i);
  auto r = std::vector<ftype>(n);
  auto rt = std::vector<ftype>(n);

  FftInternal<ftype>(a, rt, root);
  for (int i = 0; i < n; i++) r[i] = rt[i];
  FftInternal<ftype>(b, rt, root);
  for (int i = 0; i < n; i++) r[i] *= rt[i];

  // Get inverse of w
  reverse(root.begin()+1, root.end());
  FftInternal<ftype>(r, rt, root);

  std::vector<llong> res(n);
  for (int i = 0; i < n; i++) res[i] = round(rt[i].real());
  return res;
}

// FFT which produces results modulo given P.
// When "cascade" is set, the result will cascade to corresponding order, so
// that all orders greater than "cascade" will all be discarded.
std::vector<int> fft_modulo(const std::vector<int>& pa, const std::vector<int>& pb,
                            int P, int cascade = -1)
{
  int n = FindFftSize(pa, pb);
  int nP = sqrt(P)+1;
  std::vector<int> A0(pa.size(), 0), A1(pa.size(), 0);
  for (int i = 0; i < pa.size(); i++) {
    A0[i] = pa[i] % nP;
    A1[i] = pa[i] / nP;
  }
  std::vector<int> B0(pb.size(), 0), B1(pb.size(), 0);
  for (int i = 0; i < pb.size(); i++) {
    B0[i] = pb[i] % nP;
    B1[i] = pb[i] / nP;
  }

  // (A1*M+A0) * (B1*M+B0) = A1*B1*M^2 + ( (A1+A0)*(B1+B0) - A1*B1 - A0*B0 )*M + A0*B0
  //                       = A1*B1*M*(M-1) - A0*B0*(M-1) + (A1+A0)*(B1+B0)*M
  std::vector<int> res(n, 0);

  std::vector<llong> a1b1 = fft(A1, B1);
  for (int i = 0; i < n; i++) res[i] = (res[i] + a1b1[i]%P * nP%P * (nP-1)) % P;
  std::vector<llong> a0b0 = fft(A0, B0);
  for (int i = 0; i < n; i++) res[i] = (res[i] + P - a0b0[i]%P * (nP-1)%P) % P;

  for (int i = 0; i < A0.size(); i++) A0[i] += A1[i];
  for (int i = 0; i < B0.size(); i++) B0[i] += B1[i];
  std::vector<llong> asbs = fft(A0, B0);
  for (int i = 0; i < n; i++) res[i] = (res[i] + asbs[i]%P * nP) % P;

  if (cascade >= 0) {
    while (res.size() > cascade+1) res.pop_back();
  }

  return res;
}

// Fast Walsh-Hadamard transformation.
// Polynomial multiplication but with x^i (*) x^j = x^(i xor j) instead.
// Size of both a and b must be equal and are a power of 2.
template <typename T>
std::vector<T> fwht(std::vector<T> a, std::vector<T> b, FwhtOperator op)
{
  assert(a.size() == b.size());
  int n = a.size();

  FwhtInternal(a, op, false);
  FwhtInternal(b, op, false);
  std::vector<T> c(n);
  for (int i = 0; i < n; i++) c[i] = a[i]*b[i];
  FwhtInternal(c, op, true);
  return c;
}

template <typename T, typename V>
std::vector<T> fwht_pow(std::vector<T> a, V m, FwhtOperator op)
{
  int n = a.size();
  FwhtInternal(a, op, false);
  std::vector<T> c(n);
  for(int i = 0; i < n; i++) c[i] = powR(a[i], m);
  FwhtInternal(c, op, true);
  return c;
}

}  // namespae algo

#endif  // ALGO_FFT_H_
