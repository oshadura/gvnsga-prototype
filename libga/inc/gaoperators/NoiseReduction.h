#ifndef __NOISEREDUCTION__
#define __NOISEREDUCTION__

#include "generic/Population.h"

namespace geantvmoop {

template <typename Derived> class NoiseReduction {

public:
  template <typename F>
  Population<F> NoiseReduction(const Population<T> &population) {
    return static_cast<Derived *>(this)
        ->Derived::NoiseReductionImpl(population);
  }
};

} // end of namespace geantvmooop

#endif
