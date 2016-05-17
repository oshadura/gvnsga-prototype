//===--- GAValue.h - LibGA ---------------------------------*- C++
//-*-===//
//
//                     LibGA Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file GaValue.h
 * @brief Implementation of generic value class for LibGA
 * prototype
 */
//

#ifndef __GAVALUE__
#define __GAVALUE__

#include <iostream>

namespace geantvmoop {

template <typename Type> class Variable {

protected:
  Type value;

public:
  Variable() {}
  Variable(Type value) : value(value) {}

  Type getValue() const { return value; }

  virtual void setValue(const Type &value) { Variable::value = value; }
};

template <typename Type>
bool operator<(const Variable<Type> &lhs, const Variable<Type> &rhs) {
  return lhs.getValue() < rhs.getValue();
}

template <typename Type>
bool operator>(const Variable<Type> &lhs, const Variable<Type> &rhs) {
  return lhs.getValue() > rhs.getValue();
}

template <typename Type>
bool operator<=(const Variable<Type> &lhs, const Variable<Type> &rhs) {
  return lhs.getValue() <= rhs.getValue();
}

template <typename Type>
bool operator>=(const Variable<Type> &lhs, const Variable<Type> &rhs) {
  return lhs.getValue() >= rhs.getValue();
}

template <typename Type>
bool operator==(const Variable<Type> &lhs, const Variable<Type> &rhs) {
  return lhs.getValue() == rhs.getValue();
}

template <typename Type>
bool operator!=(const Variable<Type> &lhs, const Variable<Type> &rhs) {
  return lhs.getValue() != rhs.getValue();
}

template <typename Type>
std::ostream &operator<<(std::ostream &os, const Variable<Type> &rhs) {
  return os << rhs.getValue();
}
}

#endif
