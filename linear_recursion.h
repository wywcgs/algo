#ifndef ALGO_LINEAR_RECURSION_H_
#define ALGO_LINEAR_RECURSION_H_

#include <vector>

#include "defs.h"
#include "modular.h"

namespace algo {

// Berlekamp–Massey algorithm to find minimum linear recursion of given sequence g.
std::vector<int> FindMinimumLinearRecursion(const std::vector<int> g, int modp)
{
  std::vector<int> recursion, recursion_last;
  int last_value = 1, last_index = -1;

  for (int i = 0; i < g.size(); i++) {
    int v = g[i];
    for (int j = 0; j < recursion.size(); j++)
      v = (v+1LL*recursion[j]*g[i-1-j]) % modp;
    if (v == 0) continue;

    // v != 0, should fix the recursion
    llong inv = inverse(last_value, modp);
    llong c = inv*(modp-v)%modp;
    std::vector<int> t = recursion;

    int l = recursion.size();
    for (int j = 0; j <= recursion_last.size(); j++) {
      int r_index = i+j-last_index;
      while (recursion.size() < r_index)
        recursion.push_back(0);
      int term = (j == 0) ? 1 : recursion_last[j-1];
      recursion[r_index-1] += c*term%modp;
      recursion[r_index-1] %= modp;
    }

    if (2*l <= i) {
      recursion_last = t;
      last_index = i;
      last_value = v;
    }
  }

  for (int i = 0; i < recursion.size(); i++)
    recursion[i] = (modp-recursion[i]) % modp;
  reverse(recursion.begin(), recursion.end());
  return recursion;
}

// Represents n-th term with linear combination of first M terms.
// All terms with index >= M should follow the given linear recursion c
// (and c has size equal to M).
std::vector<int> FindNthElementLinearRepresenation(
    const std::vector<int>& c, llong n, int modp) {
  int m = c.size();
  if (n < m) {
    std::vector<int> r(m, 0);
    r[n] = 1;
    return r;
  }

  std::vector<int> r =
      FindNthElementLinearRepresenation(c, n/2, modp);

  std::vector<int> v(2*m, 0);

  for (int i = 0; i < m; i++) for (int j = 0; j < m; j++)
    v[i+j] = (v[i+j] + 1LL*r[i]*r[j]) % modp;
  if (n%2 == 1) {
    // Add all degress by one
    for (int i = 2*m-1; i > 0; i--) v[i] = v[i-1];
    v[0] = 0;
  }
  // Reduce higher degrees
  for (int i = 2*m-1; i >= m; i--)
    for (int j = 0; j < m; j++)
      v[i-m+j] = (v[i-m+j] + 1LL*v[i]*c[j]) % modp;

  for (int i = 0; i < m; i++) r[i] = v[i];
  return r;
}

}  // namespace algo

#endif  // ALGO_LINEAR_RECURSION_H_
