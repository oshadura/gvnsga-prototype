//===--- PolMutation.h - LibGA ---------------------------------*- C++
//-*-===//
//
//                     LibGA Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file PolMutation.h
 * @brief Implementation of  polynomial mutation for LibGA
 * prototype
 */
//===----------------------------------------------------------------------===//
#pragma once

#ifndef __GAPOLMUTATION__
#define __GAPOLMUTATION__

#include <algorithm>
#include <chrono>
#include <ctime>
#include <ostream>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "GAMutation.h"
#include "tools/Random.h"

namespace geantvmoop {

class GAPolynomialMutation : public GAMutation<GAPolynomialMutation> {

public:
  template <typename T, std::size_t SizeGene> static void MutationImpl(T &ind, double prob = -1) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> rand(0, 1);
    double fRnd, fDelta1, fDelta2, fMutPow, fDelta, fValue;
    double y, xy;
    GAVector<GADouble, SizeGene> fOffspring = static_cast<T>(ind);
    // Here we taking random 1 element..
    double lb = ind[0].GetDownLimit();
    double ub = ind[0].GetUpLimit();
    // std::cout << "Checking limits: " << lb << " and " << ub << std::endl;
    for (std::size_t j = 0; j < ind.size(); ++j) {
      if (rand(gen) <= prob) {
        y = ind[j].GetGAValue();
        fDelta1 = (y - lb) / (ub - lb);
        fDelta2 = (ub - y) / (ub - lb);
        fRnd = rand(gen);
        // CHANGE as a global parameter!
        int etamutation = 10;
        fMutPow = 1.0 / (etamutation + 1.0);
        if (fRnd <= 0.5) {
          xy = 1.0 - fDelta1;
          fValue =
              2.0 * fRnd + (1.0 - 2.0 * fRnd) * (pow(xy, (etamutation + 1.0)));
          fDelta = pow(fValue, fMutPow) - 1.0;
        } else {
          xy = 1.0 - fDelta2;
          fValue = 2.0 * (1.0 - fRnd) +
                   2.0 * (fRnd - 0.5) * (pow(xy, (etamutation + 1.0)));
          fDelta = 1.0 - (pow(fValue, fMutPow));
        }
        // std::cout << "Checking new gene: " << y << std::endl;
        y = y + fDelta * (ub - lb);
        // std::cout << "New part of Gene<T> for Mutation() " << y << std::endl;
        if (std::isinf(y)) {
          ind[j] = lb;
        }
        if (std::isnan(y)) {
          ind[j] = ub;
        }
        if (y < lb)
          y = lb;
        if (y > ub)
          y = ub;
        ind[j] = y;
      }
    }
  }
};
}

#endif
