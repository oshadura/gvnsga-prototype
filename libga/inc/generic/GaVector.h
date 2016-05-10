#ifndef __GAVECTOR__
#define __GAVECTOR__


#include <vector>

template <typename T> class GaVector : public std::vector<T> {


public:
  GaVector() : std::vector<T>() {}
  GaVector(int n, const T &v) : std::vector<T>(n, v) {}
};

#endif