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
#pragma once

#ifndef __POPULATION__
#define __POPULATION__

#include "GADouble.h"
#include "GAVector.h"
#include "TGenes.h"
#include "instrumentation/CPUManager.h"

#include <algorithm>
#include <initializer_list>
#include <iostream>
#include <list>
#include <stack>
#include <unordered_map>
#include <vector>

namespace geantvmoop {

template <typename F> class Population : public std::vector<individual_t<F> > {

public:
  // For creation of new population
  Population(std::initializer_list<individual_t<F> > list)
      : std::vector<individual_t<F> >(list) {}

  Population() : std::vector<individual_t<F> >() {}

  Population(const std::vector<individual_t<F> > &individuals)
      : std::vector<individual_t<F> >(individuals) {}

  Population(int n) {
#ifdef ENABLE_GEANTVVV
    for (int i = 0; i < n; ++i) {
      CPUManager cpumgr;
      cpumgr.InitCPU();
      hwloc_topology_t topology;
      double nbcores, ccores;
      hwloc_topology_init(&topology); // initialization
      hwloc_topology_load(topology);  // actual detection
      nbcores = hwloc_get_nbobjs_by_type(topology, HWLOC_OBJ_PU);
      hwloc_topology_destroy(topology);
      ccores =
          nbcores - cpumgr.GetCurrentValueCPU() / 100 * nbcores; // just a test
      std::cout << " Number of total free cores " << ccores << std::endl;
      if (ccores < 0.3) {
        std::cout << "Sleeping, because free CPU ratio  " << ccores
                  << " is low.." << sleep(50);
      } else {
        pid_t fArrayDead[n];
        pid_t pid = fork();
        fArrayDead[i] = pid;
        if (pid == 0) {
          typename F::Input gene = F::GetInput().random();
          auto individual = std::make_shared<TGenes<F> >(gene);
          this->push_back(individual);
          wait(NULL);
          exit(EXIT_SUCCESS);
        } else if (pid < 0) {
          std::cout << "Error on fork" << std::endl;
        }
        for (int i = 0; i < n; ++i) {
          std::cout << "Waiting for PID: " << fArrayDead[i] << " to finish.."
                    << std::endl;
          waitpid(fArrayDead[i], NULL, 0);
        }
        std::cout << "PID: " << fArrayDead[i] << " has shut down.."
                  << std::endl;
        std::fill(fArrayDead, fArrayDead + n, 0);
      }
    }
#else
    CPUManager cpumgr;
    cpumgr.InitCPU();
    hwloc_topology_t topology;
    double nbcores, ccores;
    hwloc_topology_init(&topology); // initialization
    hwloc_topology_load(topology);  // actual detection
    nbcores = hwloc_get_nbobjs_by_type(topology, HWLOC_OBJ_PU);
    hwloc_topology_destroy(topology);
    ccores =
        nbcores - cpumgr.GetCurrentValueCPU() / 100 * nbcores; // just a test
    std::cout << " Number of total free cores " << ccores << std::endl;
    if (ccores < 0.3) {
      std::cout << "Sleeping, because free CPU ratio  " << ccores << " is low.."
                << sleep(50);
    } else {
      typename F::Input gene = F::GetInput().random();
      auto individual = std::make_shared<TGenes<F> >(gene);
      this->push_back(individual);
    }
#endif
  }

  ~Population() {}

  // Stupid clang
  //#if defined __clang__
  //  void push_back(individual_t<F> ind) const { (*this).push_back(ind); }
  //#endif

  const typename F::Input &GetTGenes(int i) const {
    return (*this)[i]->GetInput();
  }

  const typename F::Output &GetTFitness(int i) const {
    return (*this)[i]->GetOutput();
  }

  typename F::Input GetGene(int gene) const {
    typename F::Input v;
    for (unsigned int j = 0; j < this->size(); ++j)
      v.push_back(GetGeneValue(j, gene));
    return v;
  }

  double GetGeneValue(int index, int gene) const {
    auto value = (*this)[index]->GetInput()[gene];
    return value.GetGAValue();
  }

  typename F::Output GetObjective(int objective) const {
    typename F::Output v;
    for (unsigned int j = 0; j < this->size(); ++j)
      v.push_back(GetObjectiveValue(j, objective));
    return v;
  }

  double GetObjectiveValue(int index, int objective) const {
    return (*this)[index]->GetOutput()[objective];
  }

  bool IsNonDominated(const individual_t<F> &ind) const {
    for (auto entry : *this) {
      if (ind->IsDominatedBy(*entry))
        return true;
    }
    return false;
  }

  const individual_t<F> PopBack() {
    individual_t<F> a = this->back();
    this->pop_back();
    return a;
  }

  void Remove(const Population<F> &pop) {
    for (auto entry : pop) {
      this->erase(std::remove(this->begin(), this->end(), entry), this->end());
    }
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
                [&m](const individual_t<F> &lhs,
                     const individual_t<F> &rhs) { return m[lhs] > m[rhs]; });
    } else
      std::sort(this->begin(), this->end(),
                [&m](const individual_t<F> &lhs,
                     const individual_t<F> &rhs) { return m[lhs] < m[rhs]; });
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
    std::cout << "---------------------------\n" << std::endl;
    for (int i = 0; i < pop.size(); ++i) {
      std::cout << "Individual " << i << std::endl;
      for (int j = 0; j < pop.GetTGenes(i).size(); ++j) {
        std::cout << pop.GetGeneValue(i, j) << "|";
      }
      std::cout << "\nFitness function value: " << std::endl;
      for (int k = 0; k < pop.GetTFitness(i).size(); ++k) {
        std::cout << pop.GetObjectiveValue(i, k) << "|";
      }
    }
    std::cout << "---------------------------\n" << std::endl;
    return s;
  }
};
}

#endif
