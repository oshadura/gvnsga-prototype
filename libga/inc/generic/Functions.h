//===--- Functions.h - LibGA ---------------------------------*- C++
//-*-===//
//
//                     LibGA Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file Functions.h
 * @brief Implementation of functions for LibGA
 * prototype
 */
//

#ifndef __FUNCTIONS__
#define __FUNCTIONS__

#include <iostream>
#include "PF.h"

namespace geantvmoop {

template <typename Derived> class Functions {

public:
  int numOfEvaluations = 0;

  static int GetNObjectives() { return Derived::GetOutput().size(); }
};
}

#endif
