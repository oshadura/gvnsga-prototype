#ifndef __SELECTION__
#define __SELECTION__

#include "generic/Population.h"

template <typename Derived> class Selection {

public:
  template <typename T>
  void SelectionGA(const Population<T> &population1, const Population<T> &population2) {
    Derived::SelectionGA(population1, population2);
  }
};

#endif
