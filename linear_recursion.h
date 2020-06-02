#ifndef ALGO_LINEAR_RECURSION_H_
#define ALGO_LINEAR_RECURSION_H_

#include <vector>

#include "defs.h"
#include "modular.h"

namespace algo {

// Berlekamp-Massey algorithm to find minimum linear recursion of given sequence g.
template<typename T>
std::vector<T> FindMinimumLinearRecursion(const std::vector<T>& g)
{
  std::vector<T> recursion, recursion_last;
  T last_value = T(1);
  int last_index = -1;

  for (int i = 0; i < g.size(); i++) {
    T v = g[i];
    for (int j = 0; j < recursion.size(); j++)
      v += recursion[j] * g[i-1-j];
    if (v == 0) continue;

    // v != 0, should fix the recursion
    std::vector<T> t = recursion;

    int l = recursion.size();
    for (int j = 0; j <= recursion_last.size(); j++) {
      int r_index = i+j-last_index;
      while (recursion.size() < r_index)
        recursion.push_back(0);
      T term = (j == 0) ? T(1) : recursion_last[j-1];
      recursion[r_index-1] -= v/last_value*term;
    }

    if (2*l <= i) {
      recursion_last = t;
      last_index = i;
      last_value = v;
    }
  }

  for (int i = 0; i < recursion.size(); i++)
    recursion[i] = -recursion[i];
  reverse(recursion.begin(), recursion.end());
  return recursion;
}

// Represents n-th term with linear combination of first M terms.
// All terms with index >= M should follow the given linear recursion c
// (and c has size equal to M).
template<typename T>
std::vector<T> FindNthElementLinearRepresenation(
    const std::vector<T>& c, llong n) {
  std::vector<bool> bin;
  for (; n != 0; n /= 2) bin.push_back(n%2 == 1);
  reverse(bin.begin(), bin.end());
  return FindNthElementLinearRepresenation(c, bin);
}

// Same with above method, but accepts n in binary form
// (higher bits first), so it can support much larger index.
template<typename T>
std::vector<T> FindNthElementLinearRepresenation(
    const std::vector<T>& c, const std::vector<bool>& bin) {
  int m = c.size();
  std::vector<T> r(m, T(0));

  int idx = 0;
  for (; idx < bin.size(); idx++)
    if (bin[idx]) break;
  if (idx == bin.size()) {
    // n = 0. result is x[0]
    r[0] = T(1);
    return r;
  }

  // starts with x[1]
  r[1] = T(1);
  for (int l = idx+1; l < bin.size(); l++) {
    std::vector<T> v(2*m, T(0));
    for (int i = 0; i < m; i++) for (int j = 0; j < m; j++)
      v[i+j] += r[i]*r[j];
    if (bin[l]) {
      // Add all degress by one
      for (int i = 2*m-1; i > 0; i--) v[i] = v[i-1];
      v[0] = T(0);
    }
    // Reduce higher degrees
    for (int i = 2*m-1; i >= m; i--) {
      if (v[i] == 0) continue;
      for (int j = 0; j < m; j++) v[i-m+j] += v[i]*c[j];
    }
    for (int i = 0; i < m; i++) r[i] = v[i];
  }

  return r;
}

}  // namespace algo

#endif  // ALGO_LINEAR_RECURSION_H_
