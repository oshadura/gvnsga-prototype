//===--- GAEvaluate.h - LibGA ---------------------------------------------*-
//C++
//-*-===//
//
//                     LibGA Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file GAEvoluate.h
 * @brief Implementation of  generic interface for evaluation for LibGA
 * prototype
 */
//===----------------------------------------------------------------------===//

#ifndef __GACROSSOVER__
#define __GACROSSOVER__

#include <memory>
#include "generic/TGenes.h"

namespace geantvmoop {

template <typename Derived> class GAEvaluate {

public:

  //void Evaluate() { output = F::Evaluate(input); }
  template <typename F> static void GAEvaluate() {
    typename F::Input individual = Derived::GAEvaluateImpl();
    return std::make_shared<geantvmoop::TGenes<F>>(individual);
  }
};
}

#endif