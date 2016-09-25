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

#include <cereal/access.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/memory.hpp>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

#include <boost/serialization/nvp.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/list.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/version.hpp>
#include <boost/serialization/optional.hpp>

namespace geantvmoop {

template <typename Type> class GAValue {

protected:
  Type value;

public:
  GAValue() {}
  GAValue(Type value) : value(value) {}
  ~GAValue() {}

  Type GetGAValue() const { return value; }

private:

  friend class cereal::access;

  //template <class Archive> void serialize(Archive &ar, TGenes<F> tg) { ar(type); }

  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive &ar, const unsigned int version) {
    ar &value;
  }

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
