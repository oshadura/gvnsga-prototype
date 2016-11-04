//===--- SBXCrossover.h - LibGA ---------------------------------*- C++
//-*-===//
//
//                     LibGA Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file GASBXCrossover.h
 * @brief Implementation of SBX crossover for LibGA
 * prototype
 */
//===----------------------------------------------------------------------===//

#ifndef __GASBXCROSSOVER__
#define __GASBXCROSSOVER__

#include "GACrossover.h"
#include "generic/GAVector.h"
#include "generic/GADouble.h"
#include "tools/Random.h"
#include <cmath>

namespace geantvmoop {

class GASBXCrossover : public GACrossover<GASBXCrossover> {

public:
  template <typename T> static T CrossoverImpl(const T &a, const T &b) {
    GAVector<GADouble> fOffspring1 = static_cast<T>(a);
    GAVector<GADouble> fOffspring2 = static_cast<T>(b);
    for (unsigned int i = 0; i < fOffspring1.size(); ++i) {
      GASBXCrossover::GeneticCrossover(fOffspring1[i], fOffspring2[i], 0.5);
    }
    return fOffspring1;
  }

  static void GeneticCrossover(GADouble &a, GADouble &b,
                               double distributionIndex) {
    double x0 = a.GetGAValue();
    double x1 = b.GetGAValue();
    double dx = fabs(x1 - x0);
    if(std::isnan(dx))
      dx = 0.000001;
    double lb = a.GetDownLimit();
    double ub = b.GetUpLimit();
    double bl;
    double bu;
    if (x0 < x1) {
      bl = 1 + 2 * (x0 - lb) / dx;
      bu = 1 + 2 * (ub - x1) / dx;
    } else {
      bl = 1 + 2 * (x1 - lb) / dx;
      bu = 1 + 2 * (ub - x0) / dx;
    }
    // use symmetric distributions
    if (bl < bu) {
      bu = bl;
    } else {
      bl = bu;
    }
    double p_bl = 1 - 1 / (2 * std::pow(bl, distributionIndex + 1));
    double p_bu = 1 - 1 / (2 * std::pow(bu, distributionIndex + 1));
    double u = Random::GetInstance().RandomDouble();
    // prevent out-of-bounds values if PRNG draws the value 1.0
    if (u == 1.0) {
      u = std::nextafter(u, -1);
    }
    double u0 = u * p_bl;
    double u1 = u * p_bu;
    double b0;
    double b1;
    if (u0 <= 0.5) {
      
      b0 = std::pow(2 * std::abs(u0), 1 / (distributionIndex + 1));
    } else {
      b0 = std::pow(0.5 / (1 - u0), 1 / (distributionIndex + 1));
    }
    if (u1 <= 0.5) {
      b1 = std::pow(2 * std::abs(u1), 1 / (distributionIndex + 1));
    } else {
      b1 = std::pow(0.5 / (1 - u1), 1 / (distributionIndex + 1));
    }
    double aValue;
    double bValue;
    if (x0 < x1) {
      aValue = 0.5 * (x0 + x1 + b0 * (x0 - x1));
      bValue = 0.5 * (x0 + x1 + b1 * (x1 - x0));
    } else {
      aValue = 0.5 * (x0 + x1 + b1 * (x0 - x1));
      bValue = 0.5 * (x0 + x1 + b0 * (x1 - x0));
    }
    // Checking nan and inf numbers
    if (std::isinf(aValue)) {
      a.SetGAValue(lb);
    }
    if (std::isnan(aValue)) {
      a.SetGAValue(ub);
    }
    if (std::isinf(bValue)) {
      b.SetGAValue(lb);
    }
    if (std::isnan(bValue)) {
      b.SetGAValue(ub);
    }
    // Guard against out-of-bounds values
    if (aValue < lb) {
      a.SetGAValue(lb);
    } else if (aValue > ub) {
      a.SetGAValue(ub);
    } else {
      a.SetGAValue(aValue);
    }
    if (bValue < lb) {
      b.SetGAValue(lb);
    } else if (bValue > ub) {
      b.SetGAValue(ub);
    } else {
      b.SetGAValue(bValue);
    }
    if (Random::GetInstance().RandomBool()) {
      double temp = a.GetGAValue();
      a.SetGAValue(b.GetGAValue());
      b.SetGAValue(temp);
    }
    // std::cout << "Crossover had been happened with " << a << " and " << b
    //          << std::endl;
  }
};
}

#endif
