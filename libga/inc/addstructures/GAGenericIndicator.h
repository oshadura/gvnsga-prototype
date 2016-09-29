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
#pragma once

#ifndef __GAGENERICINDICATOR__
#define __GAGENERICINDICATOR__

#include <algorithm>
#include <unordered_map>
#include "generic/TGenes.h"
#include "generic/Population.h"

namespace geantvmoop {

template <typename Derived, typename T> class GAGenericIndicator {

public:
  template <typename F, std::size_t SizePop>
  std::unordered_map<individual_t<F>, T> static CalculateIndicator(
      const Population<F, SizePop> &population) {
    return Derived::CalculateIndicatorImpl(population);
  }
};
}

#endif
