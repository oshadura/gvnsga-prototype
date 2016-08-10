//===--- TournamentSelection.h - LibGA ---------------------------------*- C++
//-*-===//
//
//                     LibGA Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file TournamentSelection.h
 * @brief Implementation of Tournament selection for LibGA
 * prototype
 */
//===----------------------------------------------------------------------===//
#pragma once

#ifndef __GATOURNAMENTSELECTION__
#define __GATOURNAMENTSELECTION__

#include <algorithm>
#include "GASelection.h"

namespace geantvmoop {

template <typename GAComparator>
class GATournamentSelection
    : public GASelection<GATournamentSelection<GAComparator>> {

private:
  GAComparator comp;

public:
  GATournamentSelection(const GAComparator &comp) : comp(comp) {}

  template <typename F>
  individual_t<F> UnarySelectionImpl(const Population<F> &population) {
    throw std::runtime_error("BinaryTournamentSelection does not allow to "
                             "select only single individuals!");
  }

  template <typename F>
  Population<F> MultipleSelectionImpl(const Population<F> &population, unsigned int n) {
    Population<F> result;
    Population<F> pool;
    while (result.size() < n) {
      auto index = population.GetIndex();
      std::random_shuffle(index.begin(), index.end());
      for (std::size_t i = 0; i < index.size() - 1; i += 2) {
        result.push_back(
            std::min(population[index[i]], population[index[i + 1]]));
        if (result.size() >= n)
          break;
      }
    }
    return result;
  }
};
}

#endif
