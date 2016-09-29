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
#pragma once

#ifndef __SELECTION__
#define __SELECTION__

#include "generic/Population.h"

namespace geantvmoop {

template <typename Derived> class GASelection {
public:
  template <typename F, std::size_t SizePop>
  individual_t<F> UnarySelection(const Population<F, SizePop> &population) {
    return static_cast<Derived *>(this)->UnarySelectionImpl(population);
  }

  template <typename F, std::size_t SizePop>
  Population<F, SizePop> MultipleSelection(const Population<F, SizePop> &population, int n) {
    return static_cast<Derived *>(this)->MultipleSelectionImpl(population, n);
  }
};
}

#endif
