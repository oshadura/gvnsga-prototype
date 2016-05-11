#ifndef __GACOUNTER__
#define __GACOUNTER__

#include <algorithm>
#include <unordered_map>
#include "generic/TGenes.h"
#include "generic/Population.h"

template <typename Derived, typename Type> class GACounter {

public:
  template <typename F>
  std::unordered_map<individual_t<F>, Type> static GACounterInitialize(
      const Population<F> &pop) {
    return Derived::GACounterInitialize(pop);
  }
};

#endif
