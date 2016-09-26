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
#pragma once

#ifndef __FUNCTIONS__
#define __FUNCTIONS__

#include <iostream>
#include "PF.h"

namespace geantvmoop {

template <typename Derived> class Functions {

public:
  int numOfEvaluations = 0;

private:
  static int GetNObjectives() { return Derived::GetOutput().size(); }

  template <class Archive>
  void serialize(Archive &ar, const unsigned int version) {
    ar &GetNObjectives;
  }

};
}

#endif
