#ifndef ALGO_DEFS_H_
#define ALGO_DEFS_H_

#include <functional>

namespace algo {

typedef long long llong;
typedef unsigned long long ullong;
typedef long double ldouble;

#define PB push_back
#define FI first
#define SE second
#define MP make_pair
#define CLEAR(a, v) memset((a), (v), sizeof(a))
#define ALL(v) (v).begin(), (v).end()
#define SZ(v) (int)((v).size())
#define FOR(i, a, b) for(typeof(a) i = (a); i < (b); i++)
#define FOREACH(i, v) for(const auto& i : v)

// Definitions for multiplicitive functions

typedef std::function<llong(llong, int, int)> MultiplicitiveFunction;
typedef std::function<llong(llong, int)> NtFunction;


}  // namespace algo

#endif  // ALGO_DEFS_H_
