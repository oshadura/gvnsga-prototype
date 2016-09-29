//===--- GANoiseReduction.h - LibGA ---------------------------------*- C++
//-*-===//
//
//                     LibGA Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file GANoiseReduction.h
 * @brief Implementation of noise reduction operationa for LibGA
 * prototype
 */
 //
#pragma once

#ifndef __NOISEREDUCTION__
#define __NOISEREDUCTION__

#include "generic/Population.h"

namespace geantvmoop {

template <typename Derived> class NoiseReduction {

public:
  template <typename F, std::size_t SizePop>
  Population<F, SizePop> NR(Population<F, SizePop> &population) {
    return static_cast<Derived *>(this)
        ->Derived::NoiseReductionImpl(population);
  }
};

} // end of namespace geantvmooop

#endif
