//===--- Crossover.h - LibGA --------------------------------------------*- C++
//-*-===//
//
//                     LibGA Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file Crossover.h
 * @brief Implementation of  generic interface for crossover for LibGA
 * prototype
 */
//===----------------------------------------------------------------------===//
//
#pragma once

#ifndef __GACROSSOVER__
#define __GACROSSOVER__

#include <memory>
#include "generic/TGenes.h"

namespace geantvmoop {

template <typename Derived> class GACrossover {

public:
  template <typename F>
  static individual_t<F> Crossover(individual_t<F> &a, individual_t<F> &b) {
    if (a->GetInput().size() != b->GetInput().size())
      throw std::runtime_error(
          "SinglePointCrossover is only allowed on equal sized inputs!");
    typename F::Input in = Derived::CrossoverImpl(a->GetInput(), b->GetInput());
    return std::make_shared<geantvmoop::TGenes<F>>(in);
  }
};
}

#endif
