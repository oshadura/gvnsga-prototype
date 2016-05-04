#ifndef __GACOUNTER__
#define __GACOUNTER__


#include <algorithm>
#include <unordered_map>
#include "generic/TGenes.h"
#include "generic/Population.h"

template <typename Derived, typename C> class GACounter {

public:
  template <typename F>
  std::unordered_map<std::shared_ptr<Genes<F> >, C> static GACounterInitialize(
      const Population<F> &pop) {
    return Derived::GACounterInitialize(pop);
  }
};

#endif