//===--- GACD.h - LibGA ---------------------------------*- C++
//-*-===//
//
//                     LibGA Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file GACD.h
 * @brief Implementation of crowding distance for LibGA prototype
 */
//===----------------------------------------------------------------------===//

#ifndef __GACD__
#define __GACD__

#include "GAGenericIndicator.h"

#include <vector>
#include <map>
#include <limits>
#include <cmath>
#include <algorithm>

namespace geantvmoop {

class GACD : public GAGenericIndicator<GACD, double> {

  template <typename T> using Map = std::unordered_map<individual_t<T>, double>;

public:
  template <typename T>
  static Map<T> CalculateIndicator(const Population<T> &pop) {
    std::vector<double> fMin;
    std::vector<double> fMax;
    BoundingValues(pop, fMin, fMax);
    return CalculateIndicator(pop, fMin, fMax);
  }

  template <typename T>
  static void BoundingValues(const Population<T> &pop,
                             std::vector<double> &fMin,
                             std::vector<double> &fMax) {
    for (unsigned int j = 0; j < T::GetOutput().size(); ++j) {
      auto v = pop.GetObjective(j);
      fMin.push_back(*std::min_element(v.begin(), v.end()));
      fMax.push_back(*std::max_element(v.begin(), v.end()));
    }
  }

  template <typename T>
  static Map<T> CalculateIndicator(const Population<T> &pop,
                                   std::vector<double> &fMin,
                                   std::vector<double> &fMax) {
    Map<T> fMap;
    for (auto it = pop.begin(); it != pop.end(); ++it)
      fMap[*it] = 0;
    int numOfObjectives = pop[0]->GetOutput().size();
    if (fMin.size() != numOfObjectives || fMax.size() != numOfObjectives)
      throw std::runtime_error(
          "The boundary size and objective size does not match!");
    for (int i = 0; i < numOfObjectives; ++i) {
      auto obj = pop.GetObjective(i);
      auto index = pop.SortIndex(obj);
      double fDenominator = fMax[i] - fMin[i];
      if (fDenominator < 0)
        throw std::runtime_error(
            "Error min and max values couldn't be correct!");
      fMap[pop[index[0]]] = std::numeric_limits<double>::infinity();
      fMap[pop[index[index.size() - 1]]] =
          std::numeric_limits<double>::infinity();
      for (unsigned int j = 1; j < pop.size() - 1; ++j) {
        fMap[pop[index[j]]] += (pop[index[j + 1]]->GetOutput()[i] -
                                pop[index[j - 1]]->GetOutput()[i]) /
                               fDenominator;
      }
    }
    return fMap;
  }
};
}

#endif
