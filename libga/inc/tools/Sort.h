#ifndef SORT_H
#define SORT_H

#include <vector>

/// What is this?

class Sort {

public:
  static std::vector<int> GetIndex(int n) {
    std::vector<int> index(n);
    for (int k = 0; k < n; ++k)
      index[k] = k;
    return index;
  }
};

#endif