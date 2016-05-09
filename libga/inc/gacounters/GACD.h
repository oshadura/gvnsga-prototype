#ifndef __GACD__
#define __GACD__


#include "GACounter.h"
#include <vector>
#include <map>
#include <limits>
#include <cmath>
#include <algorithm>


class GACD : public <GACD, double> {

  template <typename F>
  using Map = std::unordered_map<std::shared_ptr<Genes<F>>, double>;

public:
  template <typename F> static Map<F> CalculateCDPop(const Population<F> &pop) {
    std::vector<double> fMin;
    std::vector<double> fMax;
    CalculateBounds(pop, fMin, fMax);
    return CalculateCD(pop, fMin, fMax);
  }

  template <typename F>
  static void CalculateBounds(const Population<F> &pop, std::vector<double> &fMin,
                      std::vector<double> &fMax) {
    for (unsigned int j = 0; j < F::GetOutput().size(); ++j) {
      auto vector = pop.GetObjective(j);
      fMin.push_back(*std::min_element(vector.begin(), vector.end()));
      fMax.push_back(*std::max_element(vector.begin(), vector.end()));
    }
  }

  template <typename F>
  static Map<F> CalculateCD(const Population<F> &pop, std::vector<double> &fMin,
                           std::vector<double> &fMax) {
    typedef typename std::vector<std::shared_ptr<Genes<F>>>::iterator Iterator;
    Map<F> fMap;
    for (auto it = pop.begin(); it != pop.end(); ++it)
      fMap[*it] = 0;
    int numOfObjectives = pop[0]->getOutput().size();
    if (fMin.size() != numOfObjectives || fMax.size() != numOfObjectives)
      throw std::runtime_error(
          "The boundary size and objective size does not match!");
    for (int i = 0; i < numOfObjectives; ++i) {
      auto obj = pop.getObjective(i);
      auto index = pop.sortedIndexByVector(obj);
      double denominator = fMax[i] - fMin[i];
      if (denominator < 0)
        throw std::runtime_error(
            "Error min and max values couldn't be correct!");
      fMap[pop[index[0]]] = std::numeric_limits<double>::infinity();
      fMap[pop[index[index.size() - 1]]] = std::numeric_limits<double>::infinity();
      for (unsigned int j = 1; j < pop.size() - 1; ++j) {
        fMap[pop[index[j]]] += (pop[index[j + 1]]->getOutput()[i] -
                             pop[index[j - 1]]->getOutput()[i]) /
                            denominator;
      }
    }

    return fMap;
  }
};

#endif