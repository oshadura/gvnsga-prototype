#ifndef __GAVECTOR__
#define __GAVECTOR__

#include <vector>

namespace geantvmoop {

template <typename T> class GaVector : public std::vector<T> {

public:
  GaVector() : std::vector<T>() {}

  GaVector(int n, const T &v) : std::vector<T>(n, v) {}

  GaVector(GaVector const &);

  void operator=(GaVector const &);

  ~GaVector();

  GaVector RandomSetup() const {
    GaVector fIndividual;
    for (auto fValue : *this) {
      fIndividual.push_back(fValue.RandomSetup());
    }
    return fIndividual;
  }
};
}

#endif
