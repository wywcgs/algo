#ifndef ALGO_DEFS_H_
#define ALGO_DEFS_H_

#include <functional>
#include <utility>

namespace algo {

typedef long long llong;
typedef unsigned long long ullong;
typedef long double ldouble;
typedef std::pair<int, int> PII;

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

template <typename T>
using MultiplicitiveFunctionT = std::function<T(llong, int)>;

template <typename T>
using NtFunctionT = std::function<T(llong)>;

}  // namespace algo

#endif  // ALGO_DEFS_H_
