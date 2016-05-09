#ifndef __POPULATION__
#define __POPULATION__

#include "generic/TGenes.h"
#include "tools/Generator.h"

#include <vector>
#include <memory>
#include <algorithm>
#include <iostream>
#include <unordered_map>

template <typename F>
class Population : public std::vector<std::shared_ptr<Genes<F> > > {

  // class ContUpdatedParetoFront;

public:
  Population(std::initializer_list<std::shared_ptr<Genes<F> > > fList)
      : std::vector<std::shared_ptr<Genes<F> > >(fList) {}
  Population() : std::vector<std::shared_ptr<Genes<F> > >() {}

  Population(const std::vector<std::shared_ptr<Genes<F> > > &ind)
      : std::vector<std::shared_ptr<Genes<F> > >(ind) {}

  Population(int n) {
    // Different...
    for (int i = 0; i < n; ++i) {
      typename F::Input fRandomInd = F::GetInput().random();
      auto fIndividual = std::make_shared<Genes<F> >(fRandomInd);
      this->push_back(fIndividual);
    }
  }

  bool IsNonDominated(const std::shared_ptr<Genes<F> > &individual) const {
    for (auto fEntry : *this) {
      if (individual->IsDominated(*fEntry))
        return true;
    }
    return false;
  }

  const std::shared_ptr<Genes<F> > PopBack() {
    std::shared_ptr<Genes<F> > individual = this->back();
    this->pop_back();
    return individual;
  }

  void Remove(const Population<F> &pop) {
    for (auto fEntry : pop) {
      this->erase(std::remove(this->begin(), this->end(), fEntry), this->end());
    }
  }

  typename F::OutputType GetObjective(int fObjective) const {
    typename F::OutputType fVector;
    for (unsigned int j = 0; j < this->size(); ++j)
      fVector.push_back(GetValue(j, fObjective));
    return fVector;
  }

  double GetValue(int fIndex, int fObjective) const {
    return (*this)[fIndex]->GetOutput()[fObjective];
  }

  template <typename T>
  void SortVector(const std::vector<T> &fVector, bool IsDescending = false) {
    std::unordered_map<std::shared_ptr<Genes<F> >, T> fMap;
    for (int i = 0; i < this->size(); ++i)
      fMap[(*this)[i]] = fVector[i];
    SortMap(fMap, IsDescending);
  }

  template <typename T>
  void SortMap(std::unordered_map<std::shared_ptr<Genes<F> >, T> &fMap,
            bool IsDescending = false) {
    if (IsDescending) {
      std::sort(this->begin(), this->end(),
                [&fMap](const std::shared_ptr<Genes<F> > &lhs,
                     const std::shared_ptr<Genes<F> > &rhs) {
        return fMap[lhs] > fMap[rhs];
      });
    } else
      std::sort(this->begin(), this->end(),
                [&fMap](const std::shared_ptr<Genes<F> > &lhs,
                     const std::shared_ptr<Genes<F> > &rhs) {
        return fMap[lhs] < fMap[rhs];
      });
  }

  void SortObjective(int fObjective, bool IsDescending = false) {
    SortVector(GetObjective(fObjective), IsDescending);
  }

  template <typename T>
  std::vector<int> SortIndex(const std::vector<T> &fVector,
                             bool IsDescending = false) const {
    std::vector<int> fIndex = GetIndex();
    if (IsDescending) {
      std::sort(
          fIndex.begin(), fIndex.end(),
          [&fVector](const int &lhs, const int &rhs) { return fVector[lhs] > fVector[rhs]; });
    } else
      std::sort(
          fIndex.begin(), fIndex.end(),
          [&fVector](const int &lhs, const int &rhs) { return fVector[lhs] < fVector[rhs]; });
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
      std::cout << entry->getOutput()[0] << ", " << entry->getOutput()[1]
                << std::endl;
    }
    std::cout << "---------------------------\n";
    return s;
  }
};

#endif