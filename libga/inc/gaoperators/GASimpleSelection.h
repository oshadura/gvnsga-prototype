//===--- GASimpleSelection.h - LibGA ---------------------------------*- C++
//-*-===//
//
//                     LibGA Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file GASimpleSelection.h
 * @brief Implementation of Tournament selection for LibGA
 * prototype
 */
//===----------------------------------------------------------------------===//
#pragma once

#ifndef __GASIMPLESELECTION__
#define __GASIMPLESELECTION__

#include <algorithm>
#include "GASelection.h"
#include <random>
#include <iterator>

namespace geantvmoop {

class GASimpleSelection : public GASelection<GASimpleSelection> {

private:
  template <typename Iterator, typename RNG>
  Iterator RandomSelection(Iterator start, Iterator end, RNG &g) {
    std::uniform_int_distribution<> dis(0, std::distance(start, end) - 1);
    std::advance(start, dis(g));
    return start;
  }

  template <typename Iterator>
  Iterator RandomSalection(Iterator start, Iterator end) {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    return RandomSelection(start, end, gen);
  }

public:
  template <typename F>
  individual_t<F> UnarySelectionImpl(const Population<F> &population) {
    return *(select_randomly(population.begin(), population.end()));
  }

  template <typename F>
  Population<F> MultipleSelectionImpl(const Population<F> &population, int n) {
    Population<F> result;
    for (int i = 0; i < n; ++i) {
      result.push_back(UnarySelectionImpl(population));
    }
    return result;
  }
};
}

#endif
