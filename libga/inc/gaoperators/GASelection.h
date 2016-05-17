//===--- GASelection.h - LibGA ---------------------------------*- C++
//-*-===//
//
//                     LibGA Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file TournamentSelection.h
 * @brief Implementation of interface  for selection operator for LibGA
 * prototype
 */
//===----------------------------------------------------------------------===//

#ifndef __SELECTION__
#define __SELECTION__

#include "generic/Population.h"

namespace geantvmoop {

template <typename Derived> class GASelection {
public:
  template <typename F>
  individual_t<F> UnarySelection(const Population<F> &population) {
    return static_cast<Derived *>(this)->UnarySelectionImpl(population);
  }

  template <typename F>
  Population<F> MultipleSelection(const Population<F> &population, int n) {
    return static_cast<Derived *>(this)->MultipleSelectionImpl(population, n);
  }
};
}

#endif
