#ifndef __NOISECLEANUP__
#define __NOISECLEANUP__

#include "generic/Population.h"

namespace geantvmoop{

template <typename Derived> class NoiseCleanup {

public:
  template <typename T>
  void NoiseCleanup(const Population<T> &population1, const Population<T> &population2) {
    Derived::NoiseCleanup(population1, population2);
  }
};

}// end of namespace geantvmooop

#endif