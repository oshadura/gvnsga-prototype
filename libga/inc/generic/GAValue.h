//===--- GAValue.h - LibGA ---------------------------------*- C++
//-*-===//
//
//                     LibGA Prototype
//
//===----------------------------------------------------------------------===/

/**
 * @file GaValue.h
 * @brief Implementation of generic value class for LibGA
 * prototype
 */
 #pragma once

#ifndef __GAVALUE__
#define __GAVALUE__

#include <iostream>

namespace geantvmoop {

template <typename Type> class GAValue {

protected:
  Type value;

public:
  GAValue() {}
  GAValue(Type value) : value(value) {}
  ~GAValue() {}

  Type GetGAValue() const { return value; }

  virtual void SetGAValue(const Type &value) { GAValue::value = value; }
};

template <typename Type>
bool operator<(const GAValue<Type> &lhs, const GAValue<Type> &rhs) {
  return lhs.GetGAValue() < rhs.GetGAValue();
}

template <typename Type>
bool operator>(const GAValue<Type> &lhs, const GAValue<Type> &rhs) {
  return lhs.GetGAValue() > rhs.GetGAValue();
}

template <typename Type>
bool operator<=(const GAValue<Type> &lhs, const GAValue<Type> &rhs) {
  return lhs.GetGAValue() <= rhs.GetGAValue();
}

template <typename Type>
bool operator>=(const GAValue<Type> &lhs, const GAValue<Type> &rhs) {
  return lhs.GetGAValue() >= rhs.GetGAValue();
}

template <typename Type>
bool operator==(const GAValue<Type> &lhs, const GAValue<Type> &rhs) {
  return lhs.GetGAValue() == rhs.GetGAValue();
}

template <typename Type>
bool operator!=(const GAValue<Type> &lhs, const GAValue<Type> &rhs) {
  return lhs.GetGAValue() != rhs.GetGAValue();
}

template <typename Type>
std::ostream &operator<<(std::ostream &os, const GAValue<Type> &rhs) {
  return os << rhs.GetGAValue();
}
}

#endif
