

#ifndef MOO_SORTUTIL_H
#define MOO_SORTUTIL_H

#include <vector>

class SortUtil {

public:
  static std::vector<int> GetIndex(int n) {
    std::vector<int> index(n);
    for (int k = 0; k < n; ++k)
      index[k] = k;
    return index;
  }
};

#endif // MOO_SORTUTIL_H
