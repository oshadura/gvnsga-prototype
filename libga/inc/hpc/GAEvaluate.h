//===--- GAEvaluate.h - LibGA ---------------------------------------------*-
// C++
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
#pragma once

#ifndef __GAEVALUATE__
#define __GAEVALUATE__

#include <memory>
#include "generic/TGenes.h"

namespace geantvmoop {

template <typename Derived> class GAEvaluate {

public:
  template <typename F> static void Evaluate() { Derived::GAEvaluateImpl(); }
};
}

#endif