//===--- Population.h - LibGA ---------------------------------*- C++
//-*-===//
//
//                     LibGA Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file Population.h
 * @brief Implementation of population class for LibGA
 * prototype
 */
//

#ifndef __POPULATION__
#define __POPULATION__

#include "TGenes.h"
#include "GAVector.h"
#include "GADouble.h"
#include <vector>
#include <list>
#include <stack>
#include <iostream>
#include <algorithm>
#include <unordered_map>
#include <initializer_list>

namespace geantvmoop {

template <typename F> class Population : public std::vector<individual_t<F>> {

public:
  // Population(std::initializer_list<individual_t<F>> list)
  //    : std::vector<individual_t<F>>(list) {}

  Population() : std::vector<individual_t<F>>() {}

  Population(const std::vector<individual_t<F>> &individuals)
      : std::vector<individual_t<F>>(individuals) {}

  Population(int n) {
    for (int i = 0; i < n; ++i) {
      typename F::Input gene = F::GetInput().random();
      auto individual = std::make_shared<TGenes<F>>(gene);
      this->push_back(individual);
    }
  }

  ~Population() {}
  // Stupid clang
  void push_back(individual_t<F> ind) const { (*this).push_back(ind); }

  const typename F::Input &GetTGenes(int i) const {
    return (*this)[i]->GetInput();
  }

  const typename F::Output &GetTFitness(int i) const {
    return (*this)[i]->GetOutput();
  }

  const typename F::Input &GetIJGenes(int i, int j) const {
    return (*this)[i]->GetInput()[j];
  }

  const typename F::Output &GetIJFitness(int i, int j) const {
    return (*this)[i]->GetOutput()[j];
  }

  bool IsNonDominated(const individual_t<F> &ind) const {
    for (auto entry : *this) {
      if (ind->isDominatedBy(*entry))
        return true;
    }
    return false;
  }

  const individual_t<F> pop_back_and_delete() {
    individual_t<F> a = this->back();
    this->pop_back();
    return a;
  }

  void Remove(const Population<F> &pop) {
    for (auto entry : pop) {
      this->erase(std::remove(this->begin(), this->end(), entry), this->end());
    }
  }

  typename F::Output GetObjective(int objective) const {
    typename F::Output v;
    for (unsigned int j = 0; j < this->size(); ++j)
      v.push_back(GetGAValue(j, objective));
    return v;
  }

  double GetGAValue(int index, int objective) const {
    return (*this)[index]->GetOutput()[objective];
  }

  template <typename T>
  void SortVector(const std::vector<T> &v, bool isDescending = false) {
    std::unordered_map<individual_t<F>, T> m;
    for (int i = 0; i < this->size(); ++i)
      m[(*this)[i]] = v[i];
    SortMap(m, isDescending);
  }

  template <typename T>
  std::vector<int> SortIndex(const std::vector<T> &v,
                             bool isDescending = false) const {
    std::vector<int> index = GetIndex();
    if (isDescending) {
      std::sort(
          index.begin(), index.end(),
          [&v](const int &lhs, const int &rhs) { return v[lhs] > v[rhs]; });
    } else
      std::sort(
          index.begin(), index.end(),
          [&v](const int &lhs, const int &rhs) { return v[lhs] < v[rhs]; });
    return index;
  }

  template <typename T>
  void SortMap(std::unordered_map<individual_t<F>, T> &m,
               bool isDescending = false) {
    if (isDescending) {
      std::sort(this->begin(), this->end(),
                [&m](const individual_t<F> &lhs, const individual_t<F> &rhs) {
                  return m[lhs] > m[rhs];
                });
    } else
      std::sort(this->begin(), this->end(),
                [&m](const individual_t<F> &lhs, const individual_t<F> &rhs) {
                  return m[lhs] < m[rhs];
                });
  }

  void SortObj(int objective, bool isDescending = false) {
    SortVec(GetObjective(objective), isDescending);
  }

  std::vector<int> GetIndex() const {
    std::vector<int> index(this->size());
    for (unsigned int k = 0; k < this->size(); ++k)
      index[k] = k;
    return index;
  }

  friend std::ostream &operator<<(std::ostream &s, const Population<F> &pop) {
    std::cout << "---------------------------\n";
    std::cout << "Size of population: " << pop.size() << std::endl;
    std::cout << "---------------------------\n";
    for (int i = 0; i < pop.size(); ++i) {
      auto entry = pop[i];
      // for (int j = 0; j < entry.GetOutput().size(); ++j) {
      std::cout << entry->GetOutput()[0] << ", " << entry->GetOutput()[1]
                << ", " << entry->GetOutput()[3] << ", " << std::endl;
      // std::cout << "<- Fitness function output" << std::endl;
    }
    std::cout << "---------------------------\n";
    return s;
  }
};
}

#endif
