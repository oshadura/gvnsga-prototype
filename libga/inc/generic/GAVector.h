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
 #pragma once

#ifndef __GAVECTOR__
#define __GAVECTOR__

#include <vector>

#include <cereal/access.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/memory.hpp>

#include <boost/config.hpp>

#include <boost/archive/archive_exception.hpp>

#include <boost/archive/basic_text_iarchive.hpp>
#include <boost/archive/basic_text_oarchive.hpp>
#include <boost/archive/basic_xml_iarchive.hpp>
#include <boost/archive/basic_xml_oarchive.hpp>
#include <boost/archive/basic_binary_iarchive.hpp>
#include <boost/archive/basic_binary_oarchive.hpp>

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

template <typename Type> class GAVector : public std::vector<Type> {

public:
  GAVector() : std::vector<Type>() {}
  GAVector(int n, const Type &val) : std::vector<Type>(n, val) {}
  ~GAVector() {}

private:

  std::vector<Type> type;

  friend class cereal::access;

  //template <class Archive> void serialize(Archive &ar, TGenes<F> tg) { ar(type); }

  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive &ar, const unsigned int version) {
    ar &type;
  }

public:

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
