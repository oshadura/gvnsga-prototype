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

#ifndef __GADOUBLE__
#define __GADOUBLE__

#include "GAConstrainedValue.h"
#include "tools/Random.h"
#include <cstdlib>

namespace geantvmoop {

class GADouble : public GAConstrainedValue<double> {

public:
  GADouble(double value) : GAConstrainedValue(value, 0, 1000) {}
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
