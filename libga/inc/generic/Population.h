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
  /*
Population(std::initializer_list<std::shared_ptr<Genes<F>> > fList)
  : std::vector<std::shared_ptr<Genes<F>> >(fList) {}
  */
  Population() : std::vector<std::shared_ptr<Genes<F> > >() {}

  Population(const std::vector<std::shared_ptr<Genes<F> > > &ind)
      : std::vector<std::shared_ptr<Genes<F> > >(ind) {}

  Population(int n) {
  	//Different...
    for (int i = 0; i < n; ++i) {
      typename F::Input ind = F::GetInput().random();
      auto individual = std::make_shared<Genes<F> >(ind);
      this->push_back(individual);
    }
  }

  bool IsNonDominated(const std::shared_ptr<Genes<F> > &ind) const {
    for (auto entry : *this) {
      if (ind->IsDominated(*entry))
        return true;
    }
    return false;
  }

  const std::shared_ptr<Genes<F> > PopBack() {
    std::shared_ptr<Genes<F> > a = this->back();
    this->pop_back();
    return a;
  }

  void Remove(const Population<F> &pop) {
    for (auto entry : pop) {
      this->erase(std::remove(this->begin(), this->end(), entry), this->end());
    }
  }

  typename F::OutputType GetObjective(int objective) const {
    typename F::OutputType v;
    for (unsigned int j = 0; j < this->size(); ++j)
      v.push_back(GetValue(j, objective));
    return v;
  }

  double GetValue(int index, int objective) const {
    return (*this)[index]->GetOutput()[objective];
  }

  template <typename T>
  void Sort(const std::vector<T> &v, bool isDescending = false) {
    std::unordered_map<std::shared_ptr<Genes<F> >, T> m;
    for (int i = 0; i < this->size(); ++i)
      m[(*this)[i]] = v[i];
    Sort(m, isDescending);
  }

  template <typename T>
  void Sort(std::unordered_map<std::shared_ptr<Genes<F> >, T> &m,
            bool isDescending = false) {
    if (isDescending) {
      std::sort(this->begin(), this->end(),
                [&m](const std::shared_ptr<Genes<F> > &lhs,
                     const std::shared_ptr<Genes<F> > &rhs) {
        return m[lhs] > m[rhs];
      });
    } else
      std::sort(this->begin(), this->end(),
                [&m](const std::shared_ptr<Genes<F> > &lhs,
                     const std::shared_ptr<Genes<F> > &rhs) {
        return m[lhs] < m[rhs];
      });
  }

  void Sort(int objective, bool isDescending = false) {
    Sort(GetObjective(objective), isDescending);
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

  std::vector<int> GetIndex() const {
    std::vector<int> index(this->size());
    for (unsigned int k = 0; k < this->size(); ++k)
      index[k] = k;
    return index;
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