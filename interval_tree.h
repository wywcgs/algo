#ifndef DATA_STRUCTURE_INTERVAL_TREE_H_
#define DATA_STRUCTURE_INTERVAL_TREE_H_

#include <functional>
#include <memory>

namespace data_structure {

template <typename T>
class IntervalTree {
 public:
  IntervalTree(int x_, int y_, const T& zero_, std::function<T(T, T)> f)
      : x(x_), y(y_), data(zero_), zero(zero_), merge_f(f),
        left(nullptr), right(nullptr) {
    if (leaf()) return;
    left = std::make_unique<IntervalTree>(x, mid(), data, f);
    right = std::make_unique<IntervalTree>(mid(), y, data, f);
  }

  void Set(int px, const T& value);
  T Get(int px) const { return Get(px, px+1); }
  T Get(int start, int end) const;

 private:
  int mid() const { return (x+y)>>1; }
  bool leaf() const { return y == x+1; }

  int x;
  int y;
  T data;
  T zero;
  std::function<T(T, T)> merge_f;

  std::unique_ptr<IntervalTree<T>> left;
  std::unique_ptr<IntervalTree<T>> right;
};

template <typename T>
void IntervalTree<T>::Set(int px, const T& value) {
  if (leaf()) {
    data = value;
    return;
  }

  if (px < mid()) left->Set(px, value);
  else right->Set(px, value);

  data = merge_f(left->data, right->data);
}

template <typename T>
T IntervalTree<T>::Get(int start, int end) const {
  if (leaf()) return data;
  if (start <= x && end >= y) return data;

  T res = zero;
  int m = mid();

  if (start < m) res = merge_f(res, left->Get(start, end));
  if (end > m) res = merge_f(res, right->Get(start, end));
}

}  // namespace data_structure

#endif  // DATA_STRUCTURE_INTERVAL_TREE_H_
