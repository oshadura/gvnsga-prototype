//===--- Mutation.h - LibGA ---------------------------------*- C++
//-*-===//
//
//                     LibGA Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file Mutation.h
 * @brief Implementation of  generic interface for mutation for LibGA
 * prototype
 */
//===----------------------------------------------------------------------===//
#pragma once

#ifndef __GAMUTATION__
#define __GAMUTATION__

#include "generic/TGenes.h"

namespace geantvmoop {

template <typename Derived> class GAMutation {

public:
  template <typename F> static individual_t<F> Mutation(individual_t<F> &ind) {
    typename F::Input in = ind->GetInput();
    Derived::MutationImpl(in);
    return std::make_shared<geantvmoop::TGenes<F>>(in);
  }
};
}

#endif
