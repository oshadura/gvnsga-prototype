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

#include "PF.h"
#include <iostream>

#include <boost/archive/archive_exception.hpp>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>

#include <boost/serialization/list.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/nvp.hpp>
#include <boost/serialization/optional.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/version.hpp>

namespace geantvmoop {

template <typename Derived> class Functions {

public:
  int numOfEvaluations = 0;

private:
  static int GetNObjectives() { return Derived::GetOutput().size(); }

  friend class boost::serialization::access;

  template <class Archive>
  void serialize(Archive &ar, const unsigned int version) {
    ar &GetNObjectives;
  }
};
}

#endif
