#ifndef ALGO_NIM_H_
#define ALGO_NIM_H_

#include <algorithm>
#include <set>
#include <vector>

#include "defs.h"

namespace algo {
namespace {

std::set<llong> InitF() {
  std::set<llong> F;
  // all 2^{2^n}
  for (int i = 0; i < 6; i++)
    F.insert(1LL<<(1<<i));
  return F;
}

void FermatBase(llong a, std::vector<int>* fa) {
  for (int i = 0; (1LL<<i) <= a; i++)
    if ((a>>i)&1) fa->push_back(i);
}

}  // namespace

class NimArithmic {
 public:
  NimArithmic() : F(InitF()) {}
  llong add(llong a, llong b) const { return a^b; }
  llong multiply(llong a, llong b) const;
  llong inverse(llong a) const;

 private:
  llong MultiplyFermatBase(int na, int nb) const;

  std::set<llong> F;
};

// TODO(wywcgs): Figure out the time complexity
llong NimArithmic::multiply(llong a, llong b) const {
  if (a < b) std::swap(a, b);
  if (F.count(a) > 0) return a == b ? a/2*3 : a*b;
  if (b == 0) return 0;
  if (b == 1) return a;

  std::vector<int> fa, fb;
  FermatBase(a, &fa);
  FermatBase(b, &fb);

  llong res = 0;
  for (int na : fa) for (int nb : fb) {
    llong current = MultiplyFermatBase(na, nb);
    res = add(res, current);
  }

  return res;
}

llong NimArithmic::MultiplyFermatBase(int na, int nb) const {
  llong base = 1LL<<(na^nb);
  int overlap = na&nb;
  if (overlap == 0) return base;

  std::vector<int> bIndex;
  FermatBase(overlap, &bIndex);
  llong res = 0;
  // Define x_i := 2^{2^i}
  // x_i^2 = x_i + \prod_{j=0}^{i-1}x_j in the field
  for (int i = 0; i < (1<<bIndex.size()); i++) {
    llong current = 1;
    for (int j = 0; j < bIndex.size(); j++) {
      int thisBase = (1<<bIndex[j]) - (i>>j)%2;
      current = multiply(1<<thisBase, current);
    }
    res = add(res, current);
  }

  return multiply(res, base);
}

// Find the inverse element m of given n.
// I.e., the unique m which satisfy multiply(n, m) == 1.
// Restriction:
// - n cannot be 0, otherwise it will return -1
// - n cannot be too large (to avoid 64 bit integer overflow).
//   The exact limitation is unknown
// TODO(wywcgs): Figure out the time complexity
llong NimArithmic::inverse(llong n) const {
  if (n == 0) return -1;
  if (n == 1) return 1;

  llong gf = 2;
  for (auto f : F)
    if (f <= n) gf = f;
  llong a = n/gf, n_add_a = add(n, a);
  llong m = inverse(multiply(n, n_add_a));
  return multiply(m, n_add_a);
}

};

#endif  // ALGO_NIM_H_
