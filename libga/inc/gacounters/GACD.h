#ifndef __GACD__
#define __GACD__


#include "GACounter.h"
#include <vector>
#include <map>
#include <limits>
#include <cmath>
#include <algorithm>


class GACounter : public <CrowdingDistance, double> {

  template <typename T>
  using Map = std::unordered_map<IndividualPtr<T>, double>;

public:
  template <typename T> static Map<T> calculate_(const Population<T> &pop) {
    std::vector<double> min;
    std::vector<double> max;
    bounds_(pop, min, max);
    return calculate_(pop, min, max);
  }

  template <typename T>
  static void bounds_(const Population<T> &pop, std::vector<double> &min,
                      std::vector<double> &max) {
    for (unsigned int j = 0; j < T::getOutput().size(); ++j) {
      auto v = pop.getObjective(j);
      min.push_back(*std::min_element(v.begin(), v.end()));
      max.push_back(*std::max_element(v.begin(), v.end()));
    }
  }

  template <typename T>
  static Map<T> calculate_(const Population<T> &pop, std::vector<double> &min,
                           std::vector<double> &max) {

    typedef typename std::vector<IndividualPtr<T> >::iterator Iterator;

    Map<T> m;
    for (auto it = pop.begin(); it != pop.end(); ++it)
      m[*it] = 0;

    int numOfObjectives = pop[0]->getOutput().size();
    if (min.size() != numOfObjectives || max.size() != numOfObjectives)
      throw std::runtime_error(
          "The boundary size and objective size does not match!");

    for (int i = 0; i < numOfObjectives; ++i) {

      auto obj = pop.getObjective(i);
      auto index = pop.sortedIndexByVector(obj);

      double denominator = max[i] - min[i];
      if (denominator < 0)
        throw std::runtime_error(
            "Error min and max values couldn't be correct!");

      m[pop[index[0]]] = std::numeric_limits<double>::infinity();
      m[pop[index[index.size() - 1]]] = std::numeric_limits<double>::infinity();

      for (unsigned int j = 1; j < pop.size() - 1; ++j) {
        m[pop[index[j]]] += (pop[index[j + 1]]->getOutput()[i] -
                             pop[index[j - 1]]->getOutput()[i]) /
                            denominator;
      }
    }

    return m;
  }
};

#endif