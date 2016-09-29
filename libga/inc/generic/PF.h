//===--- PF.h - LibGA ---------------------------------*- C++
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

#ifndef __PF__
#define __PF__

#include <boost/container/static_vector.hpp>

#include "generic/Population.h"
#include <unordered_set>

namespace geantvmoop {

template <typename F, std::size_t SizePop> class PF : public std::list<individual_t<F>> {
private:
  std::list<individual_t<F>> front;

public:
  Population<F, SizePop> GetPopulation() const {
    Population<F, SizePop> result;
    for (auto ind : front)
      result.push_back(ind);
    return result;
  }

  bool Add(const individual_t<F> &ind) {
    // for every element of front
    for (auto it = front.begin(); it != front.end();) {
      // of one elements dominates ind -> does not belong to front
      if ((*it)->IsDominating(*ind) || (*it)->IsEqual(*ind))
        return false;
      // else remove all elements that are dominated by ind
      if ((*it)->IsDominated(*ind))
        front.erase(it++);
      else
        ++it;
    }
    front.push_back(ind);
    return true;
  }

  static Population<F, SizePop> ParetoFrontND(const Population<F, SizePop> &pop) {
    std::list<individual_t<F>> front;
    if (pop.empty())
      return Population<F, SizePop>();
    // function for adding an element to the front
    auto func = [&front](const individual_t<F> &ind) {
      // for every element of front
      for (auto it = front.begin(); it != front.end();) {
        // of one elements dominates ind -> does not belong to front
        if ((*it)->IsDominating(*ind))
          return false;
        // else remove all elements that are dominated by ind
        if ((*it)->IsDominated(*ind))
          front.erase(it++);
        else
          ++it;
      }
      front.push_back(ind);
      return true;
    };
    for (unsigned int i = 0; i < pop.size(); ++i)
      func(pop[i]);
    Population<F, SizePop> result;
    for (auto ind : front)
      result.push_back(ind);
    return result;
  }

  friend std::ostream &operator<<(std::ostream &os, const PF<F, SizePop> &pf) {
    for (auto i : pf)
      os << i << " ";
    std::cout << "\n";
    return os;
  }
};
}

#endif
