#ifndef __POPULATION__
#define __POPULATION__

#include "generic/TGenes.h"
#include "tools/Generator.h"

#include <vector>
#include <memory>
#include <algorithm>
#include <iostream>
#include <unordered_map>

template <typename F> class Population : public std::vector<individual_t<F>> {

  class PF;

public:
  Population(std::initializer_list<individual_t<F>> fList)
      : std::vector<individual_t<F>>(fList) {}
  Population() : std::vector<individual_t<F>>() {}

  Population(const std::vector<individual_t<F>> &ind)
      : std::vector<individual_t<F>>(ind) {}

  Population(int n) {
    // Different...
    for (int i = 0; i < n; ++i) {
      typename F::Input input = F::GetInput().RandomSetup();
      auto fInd = std::make_shared<Genes<F>>(input);
      this->push_back(fInd);
    }
  }

  bool IsNonDominated(const individual_t<F> &individual) const {
    for (auto fEntry : *this) {
      if (individual->IsDominated(*fEntry))
        return true;
    }
    return false;
  }

  const individual_t<F> PopBack() {
    individual_t<F> individual = this->back();
    this->pop_back();
    return individual;
  }

  void Remove(const Population<F> &pop) {
    for (auto fEntry : pop) {
      this->erase(std::remove(this->begin(), this->end(), fEntry), this->end());
    }
  }

  typename F::Output GetObjective(int fObjective) const {
    typename F::Output fVector;
    for (unsigned int j = 0; j < this->size(); ++j)
      fVector.push_back(GetValue(j, fObjective));
    return fVector;
  }

  double GetValue(int fIndex, int fObjective) const {
    return (*this)[fIndex]->GetOutput()[fObjective];
  }

  template <typename T>
  void SortVector(const std::vector<T> &fVector, bool IsDescending = false) {
    std::unordered_map<individual_t<F>, T> fMap;
    for (int i = 0; i < this->size(); ++i)
      fMap[(*this)[i]] = fVector[i];
    SortMap(fMap, IsDescending);
  }

  template <typename T>
  void SortMap(std::unordered_map<individual_t<F>, T> &fMap,
               bool IsDescending = false) {
    if (IsDescending) {
      std::sort(
          this->begin(), this->end(),
          [&fMap](const individual_t<F> &lhs, const individual_t<F> &rhs) {
            return fMap[lhs] > fMap[rhs];
          });
    } else
      std::sort(
          this->begin(), this->end(),
          [&fMap](const individual_t<F> &lhs, const individual_t<F> &rhs) {
            return fMap[lhs] < fMap[rhs];
          });
  }

  void SortObjective(int fObjective, bool IsDescending = false) {
    SortVector(GetObjective(fObjective), IsDescending);
  }

  template <typename T>
  std::vector<int> SortIndexVector(const std::vector<T> &fVector,
                                   bool IsDescending = false) const {
    std::vector<int> fIndex = GetIndex();
    if (IsDescending) {
      std::sort(fIndex.begin(), fIndex.end(),
                [&fVector](const int &lhs, const int &rhs) {
                  return fVector[lhs] > fVector[rhs];
                });
    } else
      std::sort(fIndex.begin(), fIndex.end(),
                [&fVector](const int &lhs, const int &rhs) {
                  return fVector[lhs] < fVector[rhs];
                });
    return fIndex;
  }

  std::vector<int> GetIndex() const {
    std::vector<int> fIndex(this->size());
    for (unsigned int k = 0; k < this->size(); ++k)
      fIndex[k] = k;
    return fIndex;
  }

  friend std::ostream &operator<<(std::ostream &s, const Population<F> &pop) {
    std::cout << "---------------------------\n";
    std::cout << "Size: " << pop.size() << std::endl;
    std::cout << "---------------------------\n";
    for (int i = 0; i < pop.size(); ++i) {
      auto entry = pop[i];
      std::cout << entry->GetOutput()[0] << ", " << entry->GetOutput()[1]
                << std::endl;
    }
    std::cout << "---------------------------\n";
    return s;
  }
};

#endif
