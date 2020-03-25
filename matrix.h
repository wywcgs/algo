#ifndef ALGO_MATRIX_H_
#define ALGO_MATRIX_H_

#include "defs.h"

namespace algo {

template<typename T, int n>
class Matrix {
 private:
  T a[n][n];

 public:
  Matrix() {
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) a[i][j] = T(0);
  }

  T* operator[](int c) { return a[c]; }
  const T* const operator[](int c) const { return a[c]; }
  Matrix<T, n>& operator+=(const Matrix<T, n>& m) { return *this = *this + m; }
  Matrix<T, n>& operator*=(const Matrix<T, n>& m) { return *this = *this * m; }
};

template<typename T, int n>
Matrix<T, n> operator+(const Matrix<T, n>& m1, const Matrix<T, n>& m2) {
  Matrix<T, n> res;
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      res[i][j] = m1[i][j] + m2[i][j];
  return std::move(res);
}

template<typename T, int n>
Matrix<T, n> operator*(const Matrix<T, n>& m1, const Matrix<T, n>& m2) {
  Matrix<T, n> res;
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      for (int k = 0; k < n; k++)
        res[i][j] += m1[i][k] * m2[k][j];
  return std::move(res);
}

// Returns m^p
template<typename T, int n>
Matrix<T, n> pow(const Matrix<T, n>& m, llong p)
{
  if (p == 0) {
    Matrix<T, n> res;
    for (int i = 0; i < n; i++) res[i][i] = T(1);
    return std::move(res);
  }

  int len = 0;
  for (llong t = p; t != 1; t /= 2) len++;

  Matrix<T, n> r = m;
  for (int i = len-1; i >= 0; i--) {
    r *= r;
    if (p&(1LL<<i)) r *= m;
  }
  return std::move(r);
}

// Return m^0 + m^1 + ... + m^p
template<typename T, int n>
Matrix<T, n> powSum(const Matrix<T, n>& m, llong p)
{
  Matrix<T, n> identity = pow(m, 0);
  if (p == 0) return identity;

  int len = 0;
  for (llong t = p; t != 1; t /= 2) len++;

  Matrix<T, n> r = m;
  for (int i = len-1; i >= 0; i--) {
    r += r*pow(m, p>>(i+1));
    if (p&(1LL<<i)) r += pow(m, p>>i);
  }
  r += identity;
  return std::move(r);
}

}  // namespace algo

#endif  // ALGO_MATRIX_H_
