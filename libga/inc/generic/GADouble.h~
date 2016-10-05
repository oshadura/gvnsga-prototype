//===--- GADouble.h - LibGA ---------------------------------*- C++
//-*-===//
//
//                     LibGA Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file GADouble.h
 * @brief Implementation of double for LibGA
 * prototype
 */
//
#pragma once

#ifndef __GADOUBLE__
#define __GADOUBLE__

#include "GAConstrainedValue.h"
#include "tools/Random.h"
#include <cstdlib>

namespace geantvmoop {

class GADouble : public GAConstrainedValue<double> {

public:
  GADouble() : GAConstrainedValue() {}
  GADouble(double value)
      : GAConstrainedValue(value, 0, 1) {} // was 0,1 // Put any other value!
  GADouble(double value, double d, double u)
      : GAConstrainedValue(value, d, u){};
  GADouble(double d, double u) : GAConstrainedValue(0, d, u){};
  ~GADouble() {}
  GADouble random() const {
    auto gene = Random::GetInstance().RandomDouble(fDown, fUp);
    return GADouble(gene, fDown, fUp);
  }
};
}

#endif
