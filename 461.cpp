#include <cstdio>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <map>
#include <set>
#include <vector>
#include <cstring>
#include <string>
#include <cmath>
#include <bitset>
#include <numeric>
#include <deque>
#include <queue>
#include <stack>
#include <memory>
#include <cassert>
#include <complex>

#include "crt.h"
#include "interval_tree.h"
#include "modular.h"
#include "roots.h"
#include "numbers.h"
#include "primes.h"
#include "division_enumerator.h"
#include "multiplicitive_sum.h"
#include "mod_num.h"
#include "defs.h"
#include "gcd.h"
#include "fft.h"
#include "linear_recursion.h"

using namespace std;
using namespace algo;
using namespace data_structure;

typedef long long ZZ;

const int N = 10000;
const ZZ E = ZZ(1e18);
const ZZ pi = E * acosl(-1);

typedef pair<ZZ, llong> pzl;

vector<ZZ> e;
vector<pzl> v;

ZZ diff = -1;
int res = 0;

bool update(ZZ d, int sum)
{
  d = abs(d);
  if (diff == -1 || diff > d) {
    diff = d;
    res = sum;
    return true;
  }
  return false;
}

void updateI(int a, int b)
{
  ZZ d = v[a].FI + v[b].FI - pi;
  int sum = v[a].SE + v[b].SE;

  update(d, sum);
}

int main()
{
  FOR(i, 0, 2*N+1) {
    ZZ k = E * (expl(1.0L*i/N) - 1);
    e.PB(k);
    if (k > pi) break;
  }

  FOR(i, 0, SZ(e)) FOR(j, i+1, SZ(e)) v.PB(pzl(e[i]+e[j], i*i+j*j));
  sort(ALL(v));

  printf("%d %d\n", SZ(e), SZ(v));

  int b = SZ(v)-1;
  FOR(i, 0, SZ(v)) {
    while (b >= 0 && v[i].FI + v[b].FI > pi) b--;
    if (b < i) break;
    updateI(i, b);
    if (b != SZ(v)-1) updateI(i, b+1);
  }

  printf("%d\n", res);

  return 0;
}

