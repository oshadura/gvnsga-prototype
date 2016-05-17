//===--- GAVector.h - LibGA ---------------------------------*- C++
//-*-===//
//
//                     LibGA Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file GaVector.h
 * @brief Implementation of vector for LibGA
 * prototype
 */
//

#ifndef __GAVECTOR__
#define __GAVECTOR__

#include <vector>

namespace geantvmoop {

template <typename Type> class GAVector : public std::vector<Type> {

public:
  GAVector() : std::vector<Type>() {}
  GAVector(int n, const Type &val) : std::vector<Type>(n, val) {}

  GAVector random() const {
    GAVector result;
    for (auto value : *this) {
      result.push_back(value.random());
    }
    return result;
  }
};
}

#endif
