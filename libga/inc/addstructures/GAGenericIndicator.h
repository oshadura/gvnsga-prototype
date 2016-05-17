//===--- GAGenericInd.h - LibGA ---------------------------------*- C++
//-*-===//
//
//                     LibGA Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file GAGenericInd.h
 * @brief Implementation of generic implementation of GA indicators for
 * algorithms for LibGA prototype
 */
//===----------------------------------------------------------------------===//

#ifndef MOO_INDICATOR_H
#define MOO_INDICATOR_H

#include <algorithm>
#include <unordered_map>
#include "generic/TGenes.h"
#include "generic/Population.h"

namespace geantvmoop {

template <typename Derived, typename T> class GAGenericIndicator {

public:
  template <typename F>
  std::unordered_map<individual_t<F>, T> static CalculateIndicator(
      const Population<F> &population) {
    return Derived::CalculateIndicator(population);
  }
};
}

#endif
