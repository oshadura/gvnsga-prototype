#ifndef __NOISECLEANUP__
#define __NOISECLEANUP__

#include "generic/Population.h"

namespace geantvmoop{

template <typename Derived> class NoiseCleanup {

public:
  template <typename F>
  Population<F> NoiseCleanup(const Population<T> &population) {
    Derived::NoiseCleanupImpl(population);
  }
};

}// end of namespace geantvmooop

#endif