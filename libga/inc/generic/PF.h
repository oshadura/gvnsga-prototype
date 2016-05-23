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

#ifndef __PF__
#define __PF__

#include "generic/Population.h"
#include <unordered_set>

namespace geantvmoop {

template <typename F> class PF {
private:
  std::list<individual_t<F>> front;

public:
  Population<F> GetPopulation() const {
    Population<F> result;
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

  static Population<F> ParetoFrontND(const Population<F> &pop) {
    std::list<individual_t<F>> front;
    if (pop.empty())
      return Population<F>();
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
    Population<F> result;
    for (auto ind : front)
      result.push_back(ind);
    return result;
  }
};
}

#endif
