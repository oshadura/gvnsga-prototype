#ifndef __SELECTION__
#define __SELECTION__

#include "generic/Population.h"

namespace geantvmoop {

template <typename Derived> class Selection {

public:
  template <typename F>
  void SelectionGABetweenPops(const Population<F> &population1,
                              const Population<F> &population2) {
    Derived::SelectionGABetweenPops(population1, population2);
  }
  template <typename F> void SelectionGAUnary(const Population<F> &population) {
    Derived::SelectionGAUnary(population);
  }
  template <typename F>
  void SelectionGAMultiple(const Population<F> &population, int n) {
    Derived::SelectionGAMultiple(population, n);
  }
};
}

#endif
