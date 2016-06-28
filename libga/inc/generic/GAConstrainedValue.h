//===--- GAConstrainedValue.h - LibGA ---------------------------------*- C++
//-*-===//
//
//                     LibGA Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file Functions.h
 * @brief Implementation of constrained value for LibGA
 * prototype
 */
//
#pragma once

#ifndef __GACONSTRAINEDVALUE__
#define __GACONSTRAINEDVALUE__

#include "GAValue.h"
#include <stdexcept>

namespace geantvmoop {

template <typename Type> class GAConstrainedValue : public GAValue<Type> {

protected:
  Type fDown;
  Type fUp;

public:
  GAConstrainedValue() {}
  GAConstrainedValue(Type d, Type u) : fDown(d), fUp(u) {};
  GAConstrainedValue(Type v, Type d, Type u)
      : GAValue<Type>(v), fDown(d), fUp(u) {};
  ~GAConstrainedValue() {}

  virtual void SetGAValue(const Type &value) {
    /*
    if (value < fDown) {
      std::cout << "Redefining value " << value << std::endl;
      this->value = fDown;
    } else if (value > fUp) {
      std::cout << "Redefining value " << value << std::endl;
      this->value = fUp;
    } else {
    */
      this->value = value;
    //}
  };

  Type GetDownLimit() const { return fDown; }

  void SetDownLimit(Type d) { fDown = d; }

  Type GetUpLimit() const { return fUp; }

  void SetUpLimit(Type u) { fUp = u; }
};
}

#endif
