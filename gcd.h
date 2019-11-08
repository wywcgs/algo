#ifndef ALGO_GCD_H_
#define ALGO_GCD_H_

namespace algo {

// For given (n1, n2, c), find a solution (x1, x2) which satisifes n1*x1 + n2*x2 = c
template <typename T>
void ExtendGcd(T n1, T n2, T c, T& x1, T& x2) {
  if (n2 == 0) {
    assert(c % n1 == 0);
  
    x1 = c / n1;
    x2 = 0;
  } else {
    T y1, y2;
    // n1*x1 + n2*x2 = c
    // ==> n2*x2 + (n1%n2 + k*n2)*x1 = c
    // ==> n2*(x2+k*x1) + n1%n2*x1 = c
    // ==> x1 = y2, x2 = y1 - k*y2
    ExtendGcd(n2, n1%n2, c, y1, y2);
    x1 = y2;
    x2 = y1 - n1/n2*y2;

    T delta = x1/n2;
    x1 -= delta*n2;
    x2 += delta*n1;
    if (x1 < 0) {
      x1 += n2;
      x2 -= n1;
    }
  }
}

}  // namespace algo

#endif  // ALGO_GCD_H_
