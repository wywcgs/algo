#ifndef ALGO_DIVISION_ENUMERATOR_H_
#define ALGO_DIVISION_ENUMERATOR_H_

#include <functional>

#include "defs.h"

namespace algo {

// Enumerates all results of pattern n/k where k iterates in [1, n].
class DivisionEnumator {
 public:
  // emit_fn(result, divisor_min, divisor_max) means:
  // - for all numbers x in range [divisor_min, divisor_max),
  //   n/x = result always holds
  //
  // It is guaranteed that all ranges forms a parition of [1, n].
  void Do(llong n, std::function<void(llong, llong, llong)> emit_fn) const;
};

void DivisionEnumator::Do(llong n, std::function<void(llong, llong, llong)> emit_fn) const {
  // Enumerates divisor
  for (llong i = 1; i*i <= n; i++) {
    llong ki = n/i;
    if (ki <= i) break;
    emit_fn(ki, i, i+1);
  }

  for (llong i = 1; i*i <= n; i++) {
    // start and end of x with floor(k/x) = i
    llong begin = n/(i+1), end = n/i;
    emit_fn(i, begin+1, end+1);
  }
}

}  // namespace algo

#endif  // ALGO_DIVISION_ENUMERATOR_H_
