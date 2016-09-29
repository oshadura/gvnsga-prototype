#pragma once

#ifndef __PCA__
#define __PCA__

#include "generic/Population.h"

namespace geantvmoop {

template <typename Derived> class PCA {
public:
  template <typename F, std::size_t SizePop> Population<F, SizePop> MVA(Population<F, SizePop> &population) {
    return static_cast<Derived *>(this)->Derived::MVAImpl(population);
  }
};
}

#endif
